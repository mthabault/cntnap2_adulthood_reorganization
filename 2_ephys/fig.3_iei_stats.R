# Load necessary libraries
library(dplyr)
library(tidyr)
library(afex)
library(emmeans)
library(readr)

# Load data
iei_data <- read.csv("C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/2_ephy/3_csv_datasets/iei.csv")

# Convert variables to factors
iei_data$genotype <- as.factor(iei_data$genotype)
iei_data$age <- as.factor(iei_data$age)
iei_data$filename <- as.factor(iei_data$filename)

# Identify bin columns
bin_columns <- names(iei_data)[!(names(iei_data) %in% c("filename", "genotype", "age"))]

# Reshape to long format
iei_long <- iei_data %>%
  pivot_longer(
    cols = all_of(bin_columns),
    names_to = "bin",
    values_to = "value"
  ) %>%
  mutate(bin = factor(bin, levels = unique(bin_columns)))

# Run RM-ANOVA function for each age
run_rm_anova <- function(data, age_label) {
  
  age_data <- data %>%
    filter(age == age_label)
  
  model <- aov_ez(
    id = "filename",
    dv = "value",
    data = age_data,
    within = "bin",
    between = "genotype",
    type = 3
  )
  
  anova_df <- as.data.frame(model$anova_table)
  anova_df$Effect <- rownames(anova_df)
  anova_df$Age <- age_label
  anova_df$Source <- "RM_ANOVA"
  
  emm <- emmeans(model, ~ genotype | bin)
  tukey_results <- contrast(emm, method = "pairwise", adjust = "tukey")
  
  tukey_df <- as.data.frame(summary(tukey_results))
  tukey_df$Age <- age_label
  tukey_df$Source <- "Tukey_posthoc"
  
  anova_out <- anova_df %>%
    mutate(
      bin = NA,
      contrast = NA
    )
  
  tukey_out <- tukey_df %>%
    mutate(
      Effect = NA
    )
  
  bind_rows(anova_out, tukey_out)
}

# Run separately for 10w and 20w
results_10w <- run_rm_anova(iei_long, "10w")
results_20w <- run_rm_anova(iei_long, "20w")

# Combine all results
combined_results <- bind_rows(results_10w, results_20w)

# Print significant RM-ANOVA effects only
sig_anova <- combined_results %>%
  filter(Source == "RM_ANOVA") %>%
  mutate(p_for_filter = coalesce(p.value, `Pr(>F)`)) %>%
  filter(!is.na(p_for_filter), p_for_filter < 0.05) %>%
  select(Age, Effect, F, p_for_filter)

cat("\nSignificant RM-ANOVA effects:\n")
print(sig_anova)

# Print significant Tukey posthoc results only
sig_tukey <- combined_results %>%
  filter(
    Source == "Tukey_posthoc",
    !is.na(p.value),
    p.value < 0.05
  ) %>%
  select(Age, bin, contrast, estimate, t.ratio, p.value)

cat("\nSignificant Tukey posthoc results:\n")
print(sig_tukey)

# Save full results to CSV
output_path <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/2_ephy/5_stat_summary/iei_stat_summary.csv"

write_csv(combined_results, output_path)

cat("\nResults saved to:\n")
cat(output_path, "\n")