library(doParallel)
library(permutes)
library(emmeans)
library(readxl)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
emm_options(pbkrtest.limit = 6000, lmerTest.limit = 6000)

data_dir      <- "/home/qmoreau/gogo_bursts/stats_new"
deriv_dir     <- "/home/common/bonaiuto/gogo_bursts/derivatives_v2/processed"
demo_dir      <- "/home/common/bonaiuto/gogo_bursts/data_v2"
clinical_xlsx <- "/home/qmoreau/gogo_bursts/stat_output/Clinical_Motor_data_for_analysis_For_Poster_With_GRVP_ZScores.xlsx"
outlier_sd    <- 3
options(contrasts = c("contr.sum", "contr.poly"))

WIN_S            <- 0.150
STEP_S           <- 0.050
MIN_CLUSTER_SIZE <- 3

# ── SHARED HELPERS ────────────────────────────────────────────────────────────

load_age_lookup <- function(path) {
  readxl::read_excel(path) %>%
    rename_with(tolower) %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    select(subject_id = study_id, age = age) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    distinct(subject_id, .keep_all = TRUE)
}

remove_rt_outliers <- function(df, sd_thresh = outlier_sd) {
  before <- nrow(df)
  df <- df %>%
    group_by(subject_id) %>%
    mutate(
      rt_mean    = mean(response_time_ms, na.rm = TRUE),
      rt_sd      = sd(response_time_ms,   na.rm = TRUE),
      is_outlier = abs(response_time_ms - rt_mean) > sd_thresh * rt_sd
    ) %>%
    ungroup() %>%
    filter(!is_outlier) %>%
    select(-rt_mean, -rt_sd, -is_outlier)
  after <- nrow(df)
  cat(sprintf("  Outliers removed: %d / %d trials (%.1f%%)\n",
              before - after, before, 100 * (before - after) / before))
  df
}

load_beta_power <- function(age_lookup, epoch) {
  library(tidyverse)
  csv_path <- file.path(data_dir, paste0("beta_power_", epoch, "_contra_trials.csv"))
  if (!file.exists(csv_path)) stop(sprintf("Run extract_beta_power_csv.py first to generate %s.", basename(csv_path)))

  df <- read_csv(csv_path, show_col_types = FALSE) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    left_join(age_lookup, by = "subject_id") %>%
    filter(!is.na(response_time), !is.na(age), response_time > 0) %>%
    mutate(
      condition        = factor(condition),
      response_time_ms = as.numeric(response_time) * 1000
    )

  time_cols <- grep("^time_", names(df), value = TRUE)
  df[time_cols] <- lapply(df[time_cols], as.numeric)
  df
}

scale_predictors <- function(df) {
  age_mean <- mean(df$age); age_sd <- sd(df$age)
  df %>% mutate(
    age_c       = (age - age_mean) / age_sd,
    trial_idx_s = scale(as.numeric(trial_idx), center = TRUE, scale = TRUE)[, 1]
  )
}

get_ratio_col <- function(df) {
  if ("z.ratio" %in% names(df)) "z.ratio" else "t.ratio"
}

build_result_row <- function(pc, factor_name, start, end, aov_tab, row_idx, em_df, ctr_df) {
  long_row  <- which(em_df$condition == "LONG")
  short_row <- which(em_df$condition == "SHORT")
  if (length(long_row) == 0 || length(short_row) == 0) return(NULL)
  ratio_col     <- get_ratio_col(em_df)
  ratio_col_ctr <- get_ratio_col(ctr_df)
  data.frame(
    factor       = factor_name,
    start        = start, end = end,
    chi_sq       = aov_tab$Chisq[row_idx],
    df           = aov_tab$Df[row_idx],
    p_val        = aov_tab[["Pr(>Chisq)"]][row_idx],
    long_z       = em_df[[ratio_col]][long_row],
    long_p       = em_df$p.value[long_row],
    short_z      = em_df[[ratio_col]][short_row],
    short_p      = em_df$p.value[short_row],
    long_short_z = ctr_df[[ratio_col_ctr]][1],
    long_short_p = ctr_df$p.value[1]
  )
}

# ── SLIDING WINDOW HELPERS ────────────────────────────────────────────────────

build_sliding_windows <- function(time_cols, win_s, step_s) {
  time_vals <- as.numeric(sub("^time_", "", time_cols))
  ord       <- order(time_vals)
  time_cols <- time_cols[ord]
  time_vals <- time_vals[ord]

  starts <- seq(min(time_vals), max(time_vals) - win_s, by = step_s)

  windows <- lapply(starts, function(s) {
    e    <- s + win_s
    cols <- time_cols[time_vals >= s & time_vals < e]
    if (length(cols) == 0) return(NULL)
    list(start = s, end = e, center = s + win_s / 2, cols = cols)
  })
  windows[!sapply(windows, is.null)]
}

compute_windowed_matrix <- function(df, windows) {
  mat <- sapply(windows, function(w) rowMeans(df[, w$cols, drop = FALSE], na.rm = TRUE))
  colnames(mat) <- paste0("w", seq_along(windows))
  mat
}

# ── WINDOW GLMM ───────────────────────────────────────────────────────────────

run_window_glmm_beta <- function(df, label, time_window) {
  library(lme4); library(car)

  window_cols <- grep("^time_", names(df), value = TRUE)
  time_vals   <- as.numeric(gsub("time_", "", window_cols))
  cols_to_avg <- window_cols[time_vals >= time_window[1] & time_vals <= time_window[2]]
  if (length(cols_to_avg) == 0) stop("No time points in the specified range.")

  df$mean_beta_power <- rowMeans(df[cols_to_avg], na.rm = TRUE)

  cat(sprintf("  Window GLMM [%s]: [%.3f, %.3f]\n", label, time_window[1], time_window[2]))

  model <- tryCatch(
    withCallingHandlers(
      glmer(
        response_time_ms ~ mean_beta_power * condition * age_c + trial_idx_s + (1 | subject_id),
        data    = df,
        family  = Gamma(link = "log"),
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      ),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) {
      cat(sprintf("  [ERROR window GLMM %s]: %s — skipping.\n", label, conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(model)) return(NULL)

  print(summary(model))
  res <- Anova(model, type = 3)
  print(res)
  pw <- emtrends(model, pairwise ~ condition, var = "mean_beta_power", infer = TRUE)
  print(pw)

  list(Anova = res, emtrends = pw)
}

# ── PERMUTATION TEST (SLIDING WINDOWS) ───────────────────────────────────────

run_trial_permutation_rt_beta <- function(epoch, win_s = WIN_S, step_s = STEP_S,
                                           min_cluster_size = MIN_CLUSTER_SIZE) {
  library(permutes); library(tidyverse)

  age_lookup      <- load_age_lookup(clinical_xlsx)
  all_win_results <- data.frame()

  df_base <- load_beta_power(age_lookup, epoch) %>%
    remove_rt_outliers() %>%
    scale_predictors()

  time_cols_all <- grep("^time_", names(df_base), value = TRUE)
  time_cols     <- time_cols_all[3:(length(time_cols_all) - 2)]

  windows <- build_sliding_windows(time_cols, win_s, step_s)
  win_starts  <- sapply(windows, function(w) w$start)
  win_ends    <- sapply(windows, function(w) w$end)
  win_centers <- sapply(windows, function(w) w$center)

  cat(sprintf("  [%s] %d sliding windows (win=%.3fs, step=%.3fs) from %d timepoints\n",
              epoch, length(windows), win_s, step_s, length(time_cols)))

  win_mat <- compute_windowed_matrix(df_base, windows)

  df_windowed <- df_base %>%
    select(subject_id, trial_idx_s, age_c, condition, response_time_ms) %>%
    bind_cols(as_tibble(win_mat))

  df_long <- df_windowed %>%
    pivot_longer(cols = starts_with("w"), names_to = "window", values_to = "beta_power") %>%
    mutate(time = win_centers[as.integer(sub("^w", "", window))]) %>%
    select(-window)

  formula <- as.formula("response_time_ms ~ beta_power * condition * age_c + trial_idx_s +
                          (1 | subject_id)")

  cl <- makeCluster(60, outfile = "")
  registerDoParallel(cl)

  perms <- clusterperm.lmer(
    formula,
    data       = df_long,
    series.var = ~time,
    nperm      = 10000,
    parallel   = TRUE,
    progress   = "",
    type       = "anova"
  )

  stopCluster(cl)
  perms$time <- as.numeric(perms$time)

  # --- Main effect of beta power ---
  main_df <- perms[perms$Factor == "beta_power" &
                     !is.na(perms$p.cluster_mass) &
                     perms$p.cluster_mass < 0.05, ]

  for (cluster_id in unique(main_df$cluster)) {
    cluster_df <- main_df[main_df$cluster == cluster_id, ]
    if (nrow(cluster_df) < min_cluster_size) next

    start <- win_starts[min(cluster_df$time)]
    end   <- win_ends[max(cluster_df$time)]

    win_results <- run_window_glmm_beta(df_base, "main", time_window = c(start, end))
    if (is.null(win_results)) next

    aov_tab  <- win_results$Anova
    main_row <- grep("^mean_beta_power$", rownames(aov_tab))
    main_p   <- aov_tab[["Pr(>Chisq)"]][main_row]

    if (main_p < 0.05) {
      em_df  <- as.data.frame(win_results$emtrends$emtrends)
      ctr_df <- as.data.frame(win_results$emtrends$contrasts)
      row <- build_result_row("beta_power", "mean_beta_power",
                              start, end, aov_tab, main_row, em_df, ctr_df)
      if (!is.null(row)) all_win_results <- rbind(all_win_results, row)
    }
  }

  # --- Interaction beta power * condition ---
  interaction_df <- perms[perms$Factor == "beta_power:condition" &
                             !is.na(perms$p.cluster_mass) &
                             perms$p.cluster_mass < 0.05, ]

  for (cluster_id in unique(interaction_df$cluster)) {
    cluster_df <- interaction_df[interaction_df$cluster == cluster_id, ]
    if (nrow(cluster_df) < min_cluster_size) next

    start <- win_starts[min(cluster_df$time)]
    end   <- win_ends[max(cluster_df$time)]

    win_results <- run_window_glmm_beta(df_base, "interaction", time_window = c(start, end))
    if (is.null(win_results)) next

    aov_tab <- win_results$Anova
    int_row <- grep("mean_beta_power:condition", rownames(aov_tab))
    int_row <- int_row[!grepl("age_c", rownames(aov_tab)[int_row])]
    if (length(int_row) == 0) next
    int_p   <- aov_tab[["Pr(>Chisq)"]][int_row]

    if (int_p < 0.05) {
      em_df  <- as.data.frame(win_results$emtrends$emtrends)
      ctr_df <- as.data.frame(win_results$emtrends$contrasts)
      row <- build_result_row("beta_power", "mean_beta_power:condition",
                              start, end, aov_tab, int_row, em_df, ctr_df)
      if (!is.null(row)) all_win_results <- rbind(all_win_results, row)
    }
  }

  write_csv(perms,           file.path(data_dir, paste0("beta_power_", epoch, "_RT_permutation.csv")))
  write_csv(all_win_results, file.path(data_dir, paste0("beta_power_", epoch, "_RT_windows.csv")))
  cat(sprintf("\nPermutation done (%s).\n", epoch))
}

# ── TREND ANALYSIS (SLIDING WINDOWS) ─────────────────────────────────────────

run_trend_analysis_rt_beta <- function(epoch, win_s = WIN_S, step_s = STEP_S) {
  library(lme4); library(lmerTest); library(tidyverse); library(emmeans)

  age_lookup    <- load_age_lookup(clinical_xlsx)
  all_trends    <- list()
  all_contrasts <- list()

  df_base <- load_beta_power(age_lookup, epoch) %>%
    remove_rt_outliers() %>%
    scale_predictors()

  time_cols_all <- grep("^time_", names(df_base), value = TRUE)
  time_cols     <- time_cols_all[3:(length(time_cols_all) - 2)]

  windows <- build_sliding_windows(time_cols, win_s, step_s)
  win_mat <- compute_windowed_matrix(df_base, windows)

  for (w_idx in seq_along(windows)) {
    time_val <- windows[[w_idx]]$center
    win_col  <- win_mat[, w_idx]
    if (sd(win_col, na.rm = TRUE) < 0.01) next
    df_col <- df_base %>% mutate(beta_power = win_col)

    model <- tryCatch(
      withCallingHandlers(
        glmer(
          response_time_ms ~ beta_power * condition * age_c + trial_idx_s + (1 | subject_id),
          data    = df_col,
          family  = Gamma(link = "log"),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        ),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) {
        cat(sprintf("  [ERROR trend %s t=%.3f]: %s — skipping.\n", epoch, time_val, conditionMessage(e)))
        return(NULL)
      }
    )

    if (is.null(model)) next

    em <- emtrends(model, pairwise ~ condition, var = "beta_power", infer = TRUE)
    tr_df <- summary(em$emtrends) %>%
      as_tibble() %>%
      mutate(time = time_val)
    ctr_df <- as.data.frame(em$contrasts) %>%
      as_tibble() %>%
      slice(1) %>%
      mutate(time = time_val)

    all_trends[[w_idx]]    <- tr_df
    all_contrasts[[w_idx]] <- ctr_df
  }

  results   <- bind_rows(all_trends) %>% rename(estimate = beta_power.trend)
  contrasts <- bind_rows(all_contrasts)

  write_csv(results,   file.path(data_dir, paste0("beta_power_", epoch, "_RT_trends_tw.csv")))
  write_csv(contrasts, file.path(data_dir, paste0("beta_power_", epoch, "_RT_contrasts_tw.csv")))
  cat(sprintf("\nTrend analysis done (%s).\n", epoch))
}

# ── RUN ───────────────────────────────────────────────────────────────────────
for (epoch in c("STIM", "RESP")) {
  run_trend_analysis_rt_beta(epoch)
  run_trial_permutation_rt_beta(epoch, min_cluster_size = MIN_CLUSTER_SIZE)
}