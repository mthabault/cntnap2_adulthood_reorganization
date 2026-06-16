library(afex)
library(tidyverse)
library(emmeans)

behaviour_data <- #paste path here

# Load data
dat <- read.csv(behaviour_data, stringsAsFactors = FALSE)

# Check factors for the tests
dat <- dat %>%
  mutate(
    mouse    = factor(mouse),
    genotype = factor(genotype),             # e.g. "wt", "ko"
    age      = factor(age, ordered = TRUE),  # e.g. "10w", "20w"
    behavior = factor(behavior)
  )

# Check structure (optional)
str(dat)

# Get list of behaviors
behaviors <- levels(dat$behavior)

anova_results <- list()

# Loop over behaviors and run rmANOVA for each
for (b in behaviors) {
  cat("\n=============================\n")
  cat("Behavior:", b, "\n")
  cat("=============================\n")
  
  dat_b <- dat %>% filter(behavior == b)
  
  # Repeated-measures ANOVA with afex::aov_ez
  aov_b <- aov_ez(
    id      = "mouse",
    dv      = "n_behaviours",
    within  = "age",
    between = "genotype",
    data    = dat_b,
    type    = 3
  )
  
  print(aov_b)
  
  # Save result object
  anova_results[[b]] <- aov_b
}

# Function to extract ANOVA table and format it
extract_anova <- function(beh) {
  tab <- anova_results[[beh]]$anova_table
  df <- as.data.frame(tab)
  df$effect <- rownames(df)
  df$behavior <- beh
  rownames(df) <- NULL
  df
}

anova_summary <- do.call(rbind, lapply(names(anova_results), extract_anova))

anova_summary <- anova_summary %>%
  select(behavior, effect, everything())

# Output path
output_path <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/data/behaviour/7_data_analysis/2_comparisons/rm_anova_results.csv"

# Save CSV
write.csv(anova_summary, output_path, row.names = FALSE)

cat("\nANOVA results saved to:\n", output_path, "\n")


posthoc_results <- list()

for (b in behaviors) {
  cat("\n=============================\n")
  cat("Behavior:", b, "\n")
  cat("=============================\n")
  
  dat_b <- dat %>% filter(behavior == b)
  
  # Run ANOVA
  aov_b <- aov_ez(
    id      = "mouse",
    dv      = "n_behaviours",
    within  = "age",
    between = "genotype",
    data    = dat_b,
    type    = 3
  )
  
  print(aov_b)
  anova_results[[b]] <- aov_b
  
  # Tukey post hoc
  emms <- emmeans(aov_b, ~ genotype * age)
  
  genotype_within_age <- contrast(
    emms, 
    method = "pairwise",
    by = "age",
    adjust = "tukey"
  ) %>%
    as.data.frame() %>%
    mutate(type = "genotype_within_age")
  
  age_within_genotype <- contrast(
    emms,
    method = "pairwise",
    by = "genotype",
    adjust = "tukey"
  ) %>%
    as.data.frame() %>%
    mutate(type = "age_within_genotype")
  
  # Combine both sets
  tuk_b <- bind_rows(genotype_within_age, age_within_genotype)
  
  tuk_b$behavior <- b
  posthoc_results[[b]] <- tuk_b
}

posthoc_summary <- bind_rows(posthoc_results)


# Reorder columns
posthoc_summary <- posthoc_summary %>%
  select(behavior, type, genotype, age, contrast, estimate, SE, df, t.ratio, p.value)

# Add significance stars
posthoc_summary <- posthoc_summary %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

print(posthoc_summary)

posthoc_path <- #paste path here

write.csv(posthoc_summary, posthoc_path, row.names = FALSE)

cat("\nTukey post-hoc results saved to:\n", posthoc_path, "\n")
