library(lme4)
library(mediation)
library(patchwork)
library(readxl)
library(tidyverse)
library(car)

data_dir      <- "/home/qmoreau/gogo_bursts/stats_new"
clinical_xlsx <- "/home/qmoreau/gogo_bursts/stat_output/Clinical_Motor_data_for_analysis_For_Poster_With_GRVP_ZScores.xlsx"

pc    <- 10
epoch <- "STIM"

# --- age lookup, same as the beta power mediation script ---
load_age_lookup <- function(path) {
  readxl::read_excel(path) %>%
    rename_with(tolower) %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    select(subject_id = study_id, age = age) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    distinct(subject_id, .keep_all = TRUE)
}

age_lookup <- load_age_lookup(clinical_xlsx)

outlier_sd <- 3

remove_rt_outliers <- function(df, sd_thresh = outlier_sd) {
  before <- nrow(df)
  df <- df %>%
    group_by(subject_id) %>%
    mutate(
      rt_mean    = mean(response_time, na.rm = TRUE),
      rt_sd      = sd(response_time,   na.rm = TRUE),
      is_outlier = abs(response_time - rt_mean) > sd_thresh * rt_sd
    ) %>%
    ungroup() %>%
    filter(!is_outlier) %>%
    select(-rt_mean, -rt_sd, -is_outlier)
  after <- nrow(df)
  cat(sprintf("  Outliers removed: %d / %d trials (%.1f%%)\n",
              before - after, before, 100 * (before - after) / before))
  df
}

run_mediation <- function(df_cond) {
  model_m <- lmer(mean_burst_count_z ~ age + (1 | subject_id), data = df_cond)
  model_y <- lmer(log(response_time) ~ mean_burst_count_z + age + (1 | subject_id), data = df_cond)
  med <- mediate(model_m, model_y, treat = "age", mediator = "mean_burst_count_z",
                 sims = 10000)
  list(med = med, model_m = model_m, model_y = model_y, s = summary(med))
}

make_arrows <- function(res, cond_label) {
  s      <- res$s
  coef_m <- summary(res$model_m)$coefficients
  coef_y <- summary(res$model_y)$coefficients

  a_est <- coef_m["age",               "Estimate"]
  a_p   <- 2 * pnorm(abs(coef_m["age",               "t value"]), lower.tail = FALSE)
  b_est <- coef_y["mean_burst_count_z", "Estimate"]
  b_p   <- 2 * pnorm(abs(coef_y["mean_burst_count_z", "t value"]), lower.tail = FALSE)

  data.frame(
    x    = c(0, 2, 0),
    y    = c(0, 1, 0),
    xend = c(2, 4, 4),
    yend = c(1, 0, 0),
    label = c(
      sprintf("a = %.2e\np = %.3f", a_est, a_p),
      sprintf("b = %.2e\np = %.3f", b_est, b_p),
      sprintf("c' = %.2e\np = %.3f", s$z0, s$z0.p)
    ),
    sig       = c(a_p < .05, b_p < .05, s$z0.p < .05),
    nudge     = c(0.08, 0.18, -0.12),
    condition = cond_label
  )
}

make_acme_box <- function(res, cond_label) {
  s <- res$s
  data.frame(
    x         = 4.3,
    y         = 1.2,
    label     = sprintf("ACME = %.2e\np = %.3f", s$d0, s$d0.p),
    sig       = s$d0.p < .05,
    condition = cond_label
  )
}

run_burst_mediation <- function(tert_val, time_window) {
  cat(sprintf("\n\n%s\nPC%d T%d %s SHORT | window: %.3f to %.3f | %s\n%s\n",
              strrep("#", 70), pc, tert_val, epoch, time_window[1], time_window[2],
              format(Sys.time()), strrep("#", 70)))

  csv_path <- file.path(data_dir, paste0("PC_", pc, "_", epoch, "_trial_burst_counts.csv"))
  if (!file.exists(csv_path)) stop(sprintf("Missing %s.", basename(csv_path)))

  df_full <- read_csv(csv_path, show_col_types = FALSE) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    left_join(age_lookup, by = "subject_id") %>%
    filter(!is.na(response_time), !is.na(age), response_time > 0) %>%
    mutate(condition = factor(condition)) %>%
    filter(tertile == tert_val, condition == "SHORT")

  window_cols <- grep("^time_", names(df_full), value = TRUE)
  time_vals   <- as.numeric(gsub("time_", "", window_cols))
  cols_to_avg <- window_cols[time_vals >= time_window[1] & time_vals <= time_window[2]]

  print(time_window)
  print(cols_to_avg)

  df_full$mean_burst_count   <- rowMeans(df_full[cols_to_avg], na.rm = TRUE)
  df_full$mean_burst_count_z <- as.numeric(scale(df_full$mean_burst_count))

  df_full <- remove_rt_outliers(df_full)

  cat(sprintf("  N trials = %d, N subjects = %d\n", nrow(df_full), n_distinct(df_full$subject_id)))

  nodes <- data.frame(
    x     = c(0, 2, 4),
    y     = c(0, 1, 0),
    label = c("Age",
              sprintf("PC%d T%d burst count\n(%.2fs to %.2fs)",
                      pc, tert_val, time_window[1], time_window[2]),
              "Response\ntime (log)")
  )

  res_short <- run_mediation(df_full)

  cat("SHORT — a:", summary(res_short$model_m)$coefficients["age","Estimate"],
      "b:", summary(res_short$model_y)$coefficients["mean_burst_count_z","Estimate"],
      "ACME:", res_short$s$d0, "p:", res_short$s$d0.p, "\n")

  arrows_short <- make_arrows(res_short, "SHORT")
  acme_short   <- make_acme_box(res_short, "SHORT")

  acme_p       <- res_short$s$d0.p
  title_str    <- sprintf("Burst count mediation - SHORT (PC%d T%d, %s)\n[%.3fs to %.3fs]",
                           pc, tert_val, epoch, time_window[1], time_window[2])
  subtitle_str <- sprintf("ACME = %.2e  |  p = %.3f", res_short$s$d0, acme_p)
  subtitle_col <- if (acme_p < .05) "#cc0000" else "grey40"

  final_plot <- ggplot() +
    geom_segment(data = arrows_short,
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = sig),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                 linewidth = 0.8) +
    geom_label(data = arrows_short,
               aes(x = (x + xend) / 2, y = (y + yend) / 2 + nudge,
                   label = label, fill = sig),
               size = 3, label.size = 0.5,
               label.padding = unit(0.3, "lines")) +
    geom_label(data = acme_short,
               aes(x = x, y = y, label = label, fill = sig),
               size = 3, label.size = 0.5, hjust = 1,
               label.padding = unit(0.3, "lines")) +
    geom_label(data = nodes,
               aes(x = x, y = y, label = label),
               size = 3, fontface = "bold", fill = "white",
               label.size = 0.5, label.padding = unit(0.4, "lines")) +
    scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid"), guide = "none") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#ffcccc"), guide = "none") +
    labs(
      title    = title_str,
      subtitle = subtitle_str
    ) +
    xlim(-0.5, 4.5) + ylim(-0.3, 1.3) +
    theme_void() +
    theme(
      plot.title        = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 6), lineheight = 1.1),
      plot.subtitle     = element_text(size = 10, hjust = 0.5, color = subtitle_col,
                                       margin = margin(b = 8),
                                       face = if (acme_p < .05) "bold" else "plain"),
      plot.background   = element_rect(fill = "white", color = NA, linewidth = 0),
      plot.margin       = margin(10, 10, 10, 10)
    )

  out_pdf <- file.path(data_dir, sprintf("burst_count_z_3SDs_PC%d_T%d_%s_SHORT_mediation_%.3f-%.3f.pdf",
                                          pc, tert_val, epoch, time_window[1], time_window[2]))
  ggsave(out_pdf, final_plot, width = 5, height = 4)
  cat(sprintf("Saved -> %s\n", out_pdf))

  cat(sprintf("\nSHORT mediation summary (PC%d T%d):\n  a (age -> burst count) = %.4e, p = %.4f\n  b (burst count -> RT)  = %.4e, p = %.4f\n  ACME                   = %.4e, p = %.4f\n",
              pc, tert_val,
              summary(res_short$model_m)$coefficients["age","Estimate"],
              2 * pnorm(abs(summary(res_short$model_m)$coefficients["age","t value"]), lower.tail = FALSE),
              summary(res_short$model_y)$coefficients["mean_burst_count_z","Estimate"],
              2 * pnorm(abs(summary(res_short$model_y)$coefficients["mean_burst_count_z","t value"]), lower.tail = FALSE),
              res_short$s$d0, res_short$s$d0.p))

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("SIMPLE LMM (SHORT, PC%d T%d): log(response_time) ~ mean_burst_count_z + (1 | subject_id)\n", pc, tert_val))
  cat(strrep("=", 70), "\n\n")

  model_short_only <- lmer(log(response_time) ~ mean_burst_count_z + (1 | subject_id),
                            data = df_full)

  cat("Type III ANOVA:\n\n")
  print(Anova(model_short_only, type = "III"))

  invisible(list(res = res_short, df = df_full, model_short_only = model_short_only))
}

# ── PC8 Tertile 1, -0.200 to 0.150 s ───────────────────────────────────────────
result_t1 <- run_burst_mediation(tert_val = 1, time_window = c(-0.200, 0.1))

# ── PC8 Tertile 3, 0.250 to 0.600 s ────────────────────────────────────────────
result_t3 <- run_burst_mediation(tert_val = 3, time_window = c(0, 0.250))