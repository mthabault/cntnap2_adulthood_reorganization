library(readr)
library(dplyr)
library(emmeans)

# Load data
data <- read.csv(#path)

# Create the group column
data$group <- paste(data$genotype, data$age, sep = "_")

# Define variables to test
variables <- c("cm", "rm", "tau")

# Prepare a dataframe to store all results
all_results <- data.frame()

# Function to perform ANOVA + Tukey and collect results
perform_anova_tukey <- function(var) {
  formula <- as.formula(paste(var, "~ group"))
  aov_model <- aov(formula, data = data)
  anova_summary <- summary(aov_model)[[1]]
  
  # Extract ANOVA stats
  anova_result <- data.frame(
    variable = var,
    test = "ANOVA",
    comparison = NA,
    Df = anova_summary["group", "Df"],
    F_value = anova_summary["group", "F value"],
    p_value = anova_summary["group", "Pr(>F)"],
    significance = ifelse(anova_summary["group", "Pr(>F)"] < 0.001, "***",
                          ifelse(anova_summary["group", "Pr(>F)"] < 0.01, "**",
                                 ifelse(anova_summary["group", "Pr(>F)"] < 0.05, "*", "ns")))
  )
  
  # Get Tukey post hoc results
  emmeans_result <- emmeans(aov_model, ~ group)
  tukey_result <- contrast(emmeans_result, method = "pairwise", adjust = "tukey")
  tukey_df <- as.data.frame(tukey_result)
  
  tukey_df$variable <- var
  tukey_df$test <- "Tukey HSD"
  tukey_df$significance <- ifelse(tukey_df$p.value < 0.001, "***",
                                  ifelse(tukey_df$p.value < 0.01, "**",
                                         ifelse(tukey_df$p.value < 0.05, "*", "ns")))
  
  # Prepare Tukey result dataframe
  tukey_result_df <- tukey_df %>%
    select(variable, test, contrast, df, estimate, SE, t.ratio, p.value, significance) %>%
    rename(comparison = contrast, Df = df, F_value = t.ratio, p_value = p.value)
  
  # Combine ANOVA and Tukey
  result_df <- bind_rows(anova_result, tukey_result_df)
  
  return(result_df)
}

# Run and collect for each variable
for (var in variables) {
  result_df <- perform_anova_tukey(var)
  all_results <- bind_rows(all_results, result_df)
}

print(all_results)

# Save all statistical results
output_path <- #path

write_csv(all_results, output_path)

cat("All statistical results saved to:", output_path, "\n")