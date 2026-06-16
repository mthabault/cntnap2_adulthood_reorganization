# Load necessary libraries
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(readr)

# Load your data
ipsc_data <- read.csv(#path)

# Genotype and age as factors
ipsc_data$genotype <- as.factor(ipsc_data$genotype)
ipsc_data$age <- as.factor(ipsc_data$age)

# Identify bin columns (= all columns are starting with X)
bin_columns <- grep("^X", names(ipsc_data), value = TRUE)

# Reshape to long format
ipsc_long <- ipsc_data %>%
  pivot_longer(cols = all_of(bin_columns),
               names_to = "bin",
               values_to = "value") %>%
  mutate(bin = factor(bin, levels = bin_columns))

# Fit linear mixed-effects model
# filename (or neuron ID) is the random effect
model <- lmer(value ~ genotype * age * bin + (1 | filename), data = ipsc_long)

# Extract ANOVA table
anova_table <- anova(model)
anova_df <- as.data.frame(anova_table)
anova_df$Effect <- rownames(anova_df)
anova_df <- anova_df %>%
  select(Effect, everything())

# Post-hoc Tukey tests for group differences at each bin
emm <- emmeans(model, ~ genotype * age | bin)
tukey_results <- contrast(emm, method = "pairwise", adjust = "tukey")

# Convert Tukey results to dataframe
tukey_df <- as.data.frame(summary(tukey_results))

# Add table labels so ANOVA and post-hoc results can be combined
anova_df <- anova_df %>%
  mutate(result_type = "LMM_ANOVA")

tukey_df <- tukey_df %>%
  mutate(result_type = "Tukey_posthoc")

# Combine LMM ANOVA and Tukey post-hoc results into one table
combined_results <- bind_rows(
  anova_df,
  tukey_df
)

# Save combined results
summary_path <- #path

write_csv(combined_results, summary_path)

cat("Combined LMM and post-hoc results saved to:", summary_path, "\n")

