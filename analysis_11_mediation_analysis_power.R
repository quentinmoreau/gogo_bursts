library(lme4)
library(mediation)
library(patchwork)
library(readxl)
library(tidyverse)

data_dir      <- "/home/qmoreau/gogo_bursts/stats_new"
clinical_xlsx <- "/home/qmoreau/gogo_bursts/stat_output/Clinical_Motor_data_for_analysis_For_Poster_With_GRVP_ZScores.xlsx"

epoch       <- "STIM"
time_window <- c(0.250, 0.600)

cat("RUN AT:", format(Sys.time()), "| window:", time_window[1], "-", time_window[2], "\n")

# --- age lookup, same as run_trend_analysis_rt_beta.R ---
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

# --- load beta power trial-level CSV, same as load_beta_power() ---
csv_path <- file.path(data_dir, "beta_power_STIM_contra_trials.csv")
if (!file.exists(csv_path)) stop("Run extract_beta_power_csv.py first to generate the CSV.")

df_full <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(subject_id = as.character(subject_id)) %>%
  left_join(age_lookup, by = "subject_id") %>%
  filter(!is.na(response_time), !is.na(age), response_time > 0) %>%
  mutate(condition = factor(condition))

window_cols <- grep("^time_", names(df_full), value = TRUE)
time_vals   <- as.numeric(gsub("time_", "", window_cols))
cols_to_avg <- window_cols[time_vals >= time_window[1] & time_vals <= time_window[2]]

print(time_window)
print(cols_to_avg)

df_full$mean_beta_power   <- rowMeans(df_full[cols_to_avg], na.rm = TRUE)
df_full$mean_beta_power_z <- as.numeric(scale(df_full$mean_beta_power))

df_full <- remove_rt_outliers(df_full)

run_mediation <- function(df_cond) {
  model_m <- lmer(mean_beta_power_z ~ age + (1 | subject_id), data = df_cond)
  model_y <- lmer(log(response_time) ~ mean_beta_power_z + age + (1 | subject_id), data = df_cond)
  med <- mediate(model_m, model_y, treat = "age", mediator = "mean_beta_power_z",
                 sims = 10000)
  list(med = med, model_m = model_m, model_y = model_y, s = summary(med))
}

make_arrows <- function(res, cond_label) {
  s      <- res$s
  coef_m <- summary(res$model_m)$coefficients
  coef_y <- summary(res$model_y)$coefficients
  
  a_est <- coef_m["age",              "Estimate"]
  a_p   <- 2 * pnorm(abs(coef_m["age",              "t value"]), lower.tail = FALSE)
  b_est <- coef_y["mean_beta_power_z", "Estimate"]
  b_p   <- 2 * pnorm(abs(coef_y["mean_beta_power_z", "t value"]), lower.tail = FALSE)
  
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

nodes <- data.frame(
  x     = c(0, 2, 4),
  y     = c(0, 1, 0),
  label = c("Age",
            sprintf("Beta power\n(%.2fs-%.2fs)",
                    time_window[1], time_window[2]),
            "Response\ntime (log)")
)

res_short <- run_mediation(filter(df_full, condition == "SHORT"))
res_long  <- run_mediation(filter(df_full, condition == "LONG"))

cat("SHORT — a:", summary(res_short$model_m)$coefficients["age","Estimate"],
    "b:", summary(res_short$model_y)$coefficients["mean_beta_power_z","Estimate"],
    "ACME:", res_short$s$d0, "p:", res_short$s$d0.p, "\n")
cat("LONG  — a:", summary(res_long$model_m)$coefficients["age","Estimate"],
    "b:", summary(res_long$model_y)$coefficients["mean_beta_power_z","Estimate"],
    "ACME:", res_long$s$d0, "p:", res_long$s$d0.p, "\n")

acme_diff    <- res_short$s$d0 - res_long$s$d0
acme_diff_se <- sqrt(res_short$s$d0.ci[2]^2 + res_long$s$d0.ci[2]^2) / (2 * qnorm(0.975))
acme_diff_z  <- acme_diff / acme_diff_se
acme_diff_p  <- 2 * pnorm(abs(acme_diff_z), lower.tail = FALSE)

cat(sprintf("Beta power mediation | ACME diff = %.3e, z = %.3f, p = %.3f\n",
            acme_diff, acme_diff_z, acme_diff_p))

arrows_all <- bind_rows(
  make_arrows(res_short, "SHORT"),
  make_arrows(res_long,  "LONG")
)

acme_boxes <- bind_rows(
  make_acme_box(res_short, "SHORT"),
  make_acme_box(res_long,  "LONG")
)

title_str    <- sprintf("Beta power mediation (%s) [%.3fs-%.3fs]", epoch, time_window[1], time_window[2])
subtitle_str <- sprintf("ACME diff = %.2e  |  p = %.3f", acme_diff, acme_diff_p)
subtitle_col <- if (acme_diff_p < .05) "#cc0000" else "grey40"

final_plot <- ggplot() +
  geom_segment(data = arrows_all,
               aes(x = x, y = y, xend = xend, yend = yend, linetype = sig),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 0.8) +
  geom_label(data = arrows_all,
             aes(x = (x + xend) / 2, y = (y + yend) / 2 + nudge,
                 label = label, fill = sig),
             size = 3, label.size = 0.5,
             label.padding = unit(0.3, "lines")) +
  geom_label(data = acme_boxes,
             aes(x = x, y = y, label = label, fill = sig),
             size = 3, label.size = 0.5, hjust = 1,
             label.padding = unit(0.3, "lines")) +
  geom_label(data = nodes,
             aes(x = x, y = y, label = label),
             size = 3, fontface = "bold", fill = "white",
             label.size = 0.5, label.padding = unit(0.4, "lines")) +
  scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid"), guide = "none") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#ffcccc"), guide = "none") +
  facet_wrap(~ condition) +
  labs(
    title    = title_str,
    subtitle = subtitle_str
  ) +
  xlim(-0.5, 4.5) + ylim(-0.3, 1.3) +
  theme_void() +
  theme(
    strip.text        = element_text(size = 12, face = "bold"),
    plot.title        = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = subtitle_col,
                                     margin = margin(b = 8),
                                     face = if (acme_diff_p < .05) "bold" else "plain"),
    plot.background   = element_rect(fill = "white", color = NA, linewidth = 0),
    plot.margin       = margin(10, 10, 10, 10)
  )

out_pdf <- file.path(data_dir, sprintf("beta_power_z_3SDs_%s_mediation_%.3f-%.3f.pdf",
                                        epoch, time_window[1], time_window[2]))
ggsave(out_pdf, final_plot, width = 8, height = 4)
cat(sprintf("Saved -> %s\n", out_pdf))
# ===== SIMPLE LMM: beta_power × condition interaction =====
library(car)

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SIMPLE LMM: log(response_time) ~ mean_beta_power_z * condition + (1 | subject_id)\n")
cat(strrep("=", 70), "\n\n")

model_interaction <- lmer(log(response_time) ~ mean_beta_power_z * condition + (1 | subject_id),
                          data = df_full)

# Print ANOVA results
cat("Type III ANOVA:\n\n")
print(Anova(model_interaction, type = "III"))


# ===== SHORT-ONLY: separate plot + stats (reuses res_short from above) =====
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SHORT ONLY — mediation plot and LMM\n")
cat(strrep("=", 70), "\n\n")

arrows_short <- filter(arrows_all, condition == "SHORT")
acme_short   <- filter(acme_boxes, condition == "SHORT")

acme_p_short       <- res_short$s$d0.p
title_str_short    <- sprintf("Beta power mediation - SHORT (%s) [%.3fs-%.3fs]", epoch, time_window[1], time_window[2])
subtitle_str_short <- sprintf("ACME = %.2e  |  p = %.3f", res_short$s$d0, acme_p_short)
subtitle_col_short <- if (acme_p_short < .05) "#cc0000" else "grey40"

short_plot <- ggplot() +
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
    title    = title_str_short,
    subtitle = subtitle_str_short
  ) +
  xlim(-0.5, 4.5) + ylim(-0.3, 1.3) +
  theme_void() +
  theme(
    plot.title        = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = subtitle_col_short,
                                     margin = margin(b = 8),
                                     face = if (acme_p_short < .05) "bold" else "plain"),
    plot.background   = element_rect(fill = "white", color = NA, linewidth = 0),
    plot.margin       = margin(10, 10, 10, 10)
  )

out_pdf_short <- file.path(data_dir, sprintf("beta_power_z_3SDs_%s_SHORT_mediation_%.3f-%.3f.pdf",
                                              epoch, time_window[1], time_window[2]))
ggsave(out_pdf_short, short_plot, width = 5, height = 4)
cat(sprintf("Saved -> %s\n", out_pdf_short))

cat(sprintf("\nSHORT mediation summary:\n  a (age -> beta) = %.4e, p = %.4f\n  b (beta -> RT)  = %.4e, p = %.4f\n  ACME            = %.4e, p = %.4f\n",
            summary(res_short$model_m)$coefficients["age","Estimate"],
            2 * pnorm(abs(summary(res_short$model_m)$coefficients["age","t value"]), lower.tail = FALSE),
            summary(res_short$model_y)$coefficients["mean_beta_power_z","Estimate"],
            2 * pnorm(abs(summary(res_short$model_y)$coefficients["mean_beta_power_z","t value"]), lower.tail = FALSE),
            res_short$s$d0, res_short$s$d0.p))

model_short_only <- lmer(log(response_time) ~ mean_beta_power_z + (1 | subject_id),
                          data = filter(df_full, condition == "SHORT"))

cat("\nSHORT-only LMM (beta power main effect, no condition term): Type III ANOVA:\n\n")
print(Anova(model_short_only, type = "III"))