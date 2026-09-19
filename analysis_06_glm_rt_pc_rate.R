library(doParallel)
library(permutes)
library(emmeans)
library(readxl)
library(lme4)
library(lmerTest)
library(tidyverse)
library(car)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
emm_options(pbkrtest.limit = 6000, lmerTest.limit = 6000)
options(contrasts = c("contr.sum", "contr.poly"))

data_dir      <- "/home/qmoreau/gogo_bursts/stats_new"
clinical_xlsx <- "/home/qmoreau/gogo_bursts/stat_output/Clinical_Motor_data_for_analysis_For_Poster_With_GRVP_ZScores.xlsx"
outlier_sd    <- 3
WIN_S             <- 0.150
STEP_S            <- 0.050
MIN_CLUSTER_SIZE  <- 3

# ── HELPERS ───────────────────────────────────────────────────────────────────

load_age_lookup <- function() {
  readxl::read_excel(clinical_xlsx) %>%
    rename_with(tolower) %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    select(subject_id = study_id, age = age) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    distinct(subject_id, .keep_all = TRUE)
}

load_tertile_df <- function(pc, epoch, tert, age_lookup) {
  csv_path <- file.path(data_dir, paste0("PC_", pc, "_", epoch, "_trial_burst_counts.csv"))
  df <- read_csv(csv_path, show_col_types = FALSE) %>%
    filter(tertile == tert) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    left_join(age_lookup, by = "subject_id") %>%
    filter(!is.na(response_time), !is.na(age), response_time > 0) %>%
    mutate(condition = factor(condition), response_time_ms = response_time * 1000)

  before <- nrow(df)
  df <- df %>%
    group_by(subject_id) %>%
    mutate(rt_mean = mean(response_time_ms), rt_sd = sd(response_time_ms)) %>%
    ungroup() %>%
    filter(abs(response_time_ms - rt_mean) <= outlier_sd * rt_sd) %>%
    select(-rt_mean, -rt_sd)
  cat(sprintf("  Outliers removed: %d / %d trials (%.1f%%)\n",
              before - nrow(df), before, 100 * (before - nrow(df)) / before))

  df %>% mutate(
    age_c       = (age - mean(age)) / sd(age),
    trial_idx_s = scale(trial_idx)[, 1]
  )
}

sliding_windows <- function(time_cols, win_s = WIN_S, step_s = STEP_S) {
  time_vals <- as.numeric(sub("^time_", "", time_cols))
  ord <- order(time_vals); time_cols <- time_cols[ord]; time_vals <- time_vals[ord]
  starts <- seq(min(time_vals), max(time_vals) - win_s, by = step_s)
  lapply(starts, function(s) {
    cols <- time_cols[time_vals >= s & time_vals < s + win_s]
    if (length(cols) == 0) return(NULL)
    list(start = s, end = s + win_s, center = s + win_s / 2, cols = cols)
  }) %>% Filter(Negate(is.null), .)
}

windowed_matrix <- function(df, windows) {
  mat <- sapply(windows, function(w) rowMeans(df[, w$cols, drop = FALSE], na.rm = TRUE))
  colnames(mat) <- paste0("w", seq_along(windows))
  mat
}

fit_glmm <- function(df, label, time_val = NA) {
  tryCatch(
    withCallingHandlers(
      glmer(response_time_ms ~ burst_count * condition * age_c + trial_idx_s + (1 | subject_id),
            data = df, family = Gamma(link = "log"),
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) {
      cat(sprintf("  [ERROR %s t=%s]: %s — skipping.\n", label, time_val, conditionMessage(e)))
      NULL
    }
  )
}

ratio_col <- function(df) if ("z.ratio" %in% names(df)) "z.ratio" else "t.ratio"

result_row <- function(factor_name, start, end, aov_tab, row_idx, em_df, ctr_df, pc, tert) {
  long_row <- which(em_df$condition == "LONG"); short_row <- which(em_df$condition == "SHORT")
  if (length(long_row) == 0 || length(short_row) == 0) return(NULL)
  rc <- ratio_col(em_df); rc_ctr <- ratio_col(ctr_df)
  data.frame(
    pc = pc, tertile = tert, factor = factor_name, start = start, end = end,
    chi_sq = aov_tab$Chisq[row_idx], df = aov_tab$Df[row_idx], p_val = aov_tab[["Pr(>Chisq)"]][row_idx],
    long_z = em_df[[rc]][long_row], long_p = em_df$p.value[long_row],
    short_z = em_df[[rc]][short_row], short_p = em_df$p.value[short_row],
    long_short_z = ctr_df[[rc_ctr]][1], long_short_p = ctr_df$p.value[1]
  )
}

# ── TREND ANALYSIS (sliding window GLMM) ──────────────────────────────────────

run_trend <- function(pc, epoch, tert, df) {
  time_cols <- grep("^time_", names(df), value = TRUE)
  time_cols <- time_cols[3:(length(time_cols) - 2)]
  windows   <- sliding_windows(time_cols)
  win_mat   <- windowed_matrix(df, windows)

  trends <- list(); contrasts <- list()

  for (w in seq_along(windows)) {
    if (sd(win_mat[, w], na.rm = TRUE) < 0.01) next
    df_w <- df %>% mutate(burst_count = win_mat[, w])
    model <- fit_glmm(df_w, sprintf("PC%d T%d %s", pc, tert, epoch), windows[[w]]$center)
    if (is.null(model)) next

    em <- emtrends(model, pairwise ~ condition, var = "burst_count", infer = TRUE)
    trends[[w]]    <- summary(em$emtrends) %>% as_tibble() %>% mutate(time = windows[[w]]$center, tertile = tert, pc = pc)
    contrasts[[w]] <- as.data.frame(em$contrasts) %>% as_tibble() %>% slice(1) %>% mutate(time = windows[[w]]$center, tertile = tert, pc = pc)
  }

  if (length(trends) == 0) return(NULL)
  list(trends = bind_rows(trends) %>% rename(estimate = burst_count.trend),
       contrasts = bind_rows(contrasts))
}

# ── PERMUTATION (sliding window cluster test) ─────────────────────────────────

run_permutation <- function(pc, epoch, tert, df) {
  time_cols <- grep("^time_", names(df), value = TRUE)
  time_cols <- time_cols[3:(length(time_cols) - 2)]
  windows   <- sliding_windows(time_cols)
  win_mat   <- windowed_matrix(df, windows)
  win_starts  <- sapply(windows, `[[`, "start")
  win_ends    <- sapply(windows, `[[`, "end")
  win_centers <- sapply(windows, `[[`, "center")

  df_long <- df %>%
    select(subject_id, trial_idx_s, age_c, condition, response_time_ms) %>%
    bind_cols(as_tibble(win_mat)) %>%
    pivot_longer(starts_with("w"), names_to = "window", values_to = "burst_count") %>%
    mutate(time = win_centers[as.integer(sub("^w", "", window))]) %>%
    select(-window)

  cl <- makeCluster(60, outfile = "")
  registerDoParallel(cl)
  perms <- clusterperm.lmer(
    response_time_ms ~ burst_count * condition * age_c + trial_idx_s + (1 | subject_id),
    data = df_long, series.var = ~time, nperm = 10000, parallel = TRUE, progress = "", type = "anova"
  )
  stopCluster(cl)
  perms$time <- as.numeric(perms$time)
  perms$tertile <- tert
  perms$pc <- pc

  win_rows <- data.frame()

  for (factor_name in c("burst_count", "burst_count:condition")) {
    clu_df <- perms[perms$Factor == factor_name & !is.na(perms$p.cluster_mass) & perms$p.cluster_mass < 0.05, ]
    for (cid in unique(clu_df$cluster)) {
      cdf <- clu_df[clu_df$cluster == cid, ]
      if (nrow(cdf) < MIN_CLUSTER_SIZE) next
      start <- win_starts[min(cdf$time)]; end <- win_ends[max(cdf$time)]

      cols_to_avg <- time_cols[as.numeric(sub("^time_", "", time_cols)) >= start &
                                as.numeric(sub("^time_", "", time_cols)) <= end]
      df_win <- df %>% mutate(burst_count = rowMeans(df[cols_to_avg], na.rm = TRUE))
      model <- fit_glmm(df_win, sprintf("PC%d T%d %s [%s]", pc, tert, epoch, factor_name))
      if (is.null(model)) next

      aov_tab <- Anova(model, type = 3)
      row_idx <- grep(paste0("^", gsub(":", ":", factor_name), "$"), rownames(aov_tab))
      if (length(row_idx) == 0 || aov_tab[["Pr(>Chisq)"]][row_idx] >= 0.05) next

      em <- emtrends(model, pairwise ~ condition, var = "burst_count", infer = TRUE)
      row <- result_row(factor_name, start, end, aov_tab, row_idx,
                        as.data.frame(em$emtrends), as.data.frame(em$contrasts), pc, tert)
      if (!is.null(row)) win_rows <- rbind(win_rows, row)
    }
  }

  list(perms = perms, windows = win_rows)
}

# ── RUN ───────────────────────────────────────────────────────────────────────

age_lookup <- load_age_lookup()

for (pc in c(8, 10, 12)) {
  for (epoch in c("STIM", "RESP")) {
    trends <- list(); contrasts <- list()
    perms  <- list(); win_results <- list()

    for (tert in 1:3) {
      df <- load_tertile_df(pc, epoch, tert, age_lookup)

      tr <- run_trend(pc, epoch, tert, df)
      if (!is.null(tr)) { trends[[tert]] <- tr$trends; contrasts[[tert]] <- tr$contrasts }

      pm <- run_permutation(pc, epoch, tert, df)
      perms[[tert]] <- pm$perms
      win_results[[tert]] <- pm$windows
    }

    if (length(trends) > 0) {
      write_csv(bind_rows(trends),    file.path(data_dir, paste0("PC_", pc, "_", epoch, "_RT_trends_tw.csv")))
      write_csv(bind_rows(contrasts), file.path(data_dir, paste0("PC_", pc, "_", epoch, "_RT_contrasts_tw.csv")))
    }
    write_csv(bind_rows(perms),       file.path(data_dir, paste0("PC_", pc, "_", epoch, "_RT_permutation_tw.csv")))
    write_csv(bind_rows(win_results), file.path(data_dir, paste0("PC_", pc, "_", epoch, "_RT_windows_tw.csv")))

    cat(sprintf("Done: PC%d %s\n", pc, epoch))
  }
}