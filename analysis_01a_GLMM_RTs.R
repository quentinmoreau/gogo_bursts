library(lme4)
library(lmerTest)
library(tidyverse)
library(emmeans)
library(car)
library(readxl)

emm_options(pbkrtest.limit = 20000, lmerTest.limit = 20000)

data_dir      <- "/home/qmoreau/gogo_bursts/stats_new"
pc            <- "12"
epoch         <- "STIM"
clinical_xlsx <- "/home/qmoreau/Downloads/Clinical_Motor_data_for_analysis_For_Poster_With_GRVP_ZScores.xlsx"
outlier_sd    <- 3

tag <- sprintf("PC_%s_%s_RTage", pc, epoch)

# ── LOAD AGE ──────────────────────────────────────────────────────────────────
age_lookup <- read_excel(clinical_xlsx) %>%
  rename_with(tolower) %>%
  rename_with(~ str_replace_all(.x, " ", "_")) %>%
  select(subject_id = study_id, age = age) %>%
  mutate(subject_id = as.character(subject_id)) %>%
  distinct(subject_id, .keep_all = TRUE)

cat(sprintf("Age lookup: %d subjects, age range %.1f - %.1f\n",
            nrow(age_lookup), min(age_lookup$age), max(age_lookup$age)))

# ── LOAD TRIAL DATA ───────────────────────────────────────────────────────────
df_raw <- read_csv(paste0(data_dir, "/PC_", pc, "_", epoch, "_trial_burst_counts.csv")) %>%
  mutate(subject_id = as.character(subject_id)) %>%
  left_join(age_lookup, by = "subject_id") %>%
  filter(!is.na(response_time), !is.na(age), response_time > 0) %>%
  select(subject_id, trial_idx, condition, response_time, age) %>%
  mutate(
    condition        = factor(condition),
    response_time_ms = response_time * 1000
  )

# ── OUTLIER REMOVAL ───────────────────────────────────────────────────────────
df <- df_raw %>%
  group_by(subject_id) %>%
  mutate(
    rt_mean    = mean(response_time_ms),
    rt_sd      = sd(response_time_ms),
    is_outlier = abs(response_time_ms - rt_mean) > outlier_sd * rt_sd
  ) %>%
  ungroup()

n_outliers <- sum(df$is_outlier)
cat(sprintf("Outliers removed: %d / %d trials (%.1f%%)\n",
            n_outliers, nrow(df), 100 * n_outliers / nrow(df)))

df <- df %>%
  filter(!is_outlier) %>%
  select(-rt_mean, -rt_sd, -is_outlier)

age_mean <- mean(df$age)
age_sd   <- sd(df$age)

df <- df %>%
  mutate(
    age_c       = (age - age_mean) / age_sd,
    trial_idx_s = scale(trial_idx, center = TRUE, scale = TRUE)[, 1]
  )

cat(sprintf("df after removal: %d rows, %d subjects, age range %.1f - %.1f\n",
            nrow(df), n_distinct(df$subject_id),
            min(df$age), max(df$age)))

# ── GAMMA GLMM ───────────────────────────────────────────────────────────────
cat("\n======================================\n")
cat(sprintf("Gamma GLMM (log link) | outlier threshold: ±%.1f SD\n", outlier_sd))
cat("======================================\n")

m_gamma <- glmer(
  response_time_ms ~ condition * age_c + trial_idx_s + (1 | subject_id),
  data    = df,
  family  = Gamma(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(m_gamma))

anova_tbl <- Anova(m_gamma, type = 3)
cat("\n=== Type III Wald ===\n")
print(anova_tbl)

cat("\n=== Age trends by condition ===\n")
age_trends <- emtrends(m_gamma, pairwise ~ condition, var = "age_c", infer = TRUE)
print(age_trends)

as.data.frame(age_trends$emtrends) %>%
  mutate(slope_per_year = age_c.trend / age_sd) %>%
  write_csv(file.path(data_dir, paste0(tag, "_emtrends.csv")))

as.data.frame(age_trends$contrasts) %>%
  write_csv(file.path(data_dir, paste0(tag, "_emtrends_contrasts.csv")))

cat("\n=== EMMs (back-transformed to ms) ===\n")
emm <- emmeans(m_gamma, ~ condition, type = "response")
print(emm)
print(pairs(emm, adjust = "bonferroni"))

as.data.frame(anova_tbl) %>%
  rownames_to_column("term") %>%
  write_csv(file.path(data_dir, paste0(tag, "_anova.csv")))

# ── PREDICTIONS ───────────────────────────────────────────────────────────────
age_seq <- seq(min(df$age), max(df$age), length.out = 100)

pred_raw <- emmeans(m_gamma, ~ condition,
                    at = list(age_c = (age_seq - age_mean) / age_sd),
                    by = "age_c",
                    type = "response") %>%
  as.data.frame() %>%
  mutate(age = age_c * age_sd + age_mean)

pred_df <- pred_raw %>% select(age, condition, response, SE, asymp.LCL, asymp.UCL)

write_csv(pred_df, file.path(data_dir, paste0(tag, "_predictions.csv")))

obs_df <- df %>%
  group_by(subject_id, condition, age) %>%
  summarise(rt_ms = mean(response_time_ms, na.rm = TRUE), .groups = "drop")

as.data.frame(age_trends$emtrends) %>%
  mutate(
    slope_per_year = age_c.trend / age_sd,
    pct_per_year   = 100 * (exp(slope_per_year) - 1)
  ) %>%
  write_csv(file.path(data_dir, paste0(tag, "_emtrends.csv")))


cat(sprintf("Wrote %s_*.csv to %s\n", tag, data_dir))