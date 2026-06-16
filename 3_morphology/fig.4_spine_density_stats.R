library(readr)
library(dplyr)
library(emmeans)

spines <- read_csv("C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/3_morphology/3_csv_datasets/spines_mastersheet.csv")

dendrites <- read_csv("C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/3_morphology/3_csv_datasets/spines_dendrites_info.csv")

dendrites_len <- dendrites %>%
  select(mouse_id, cell_id, ID, dendrite_length) %>%
  rename(
    parent_id = ID,
    segment_length = dendrite_length
  )

spines_seg <- spines %>%
  filter(depth %in% c(1, 2, 3, 4)) %>%
  left_join(
    dendrites_len,
    by = c("mouse_id", "cell_id", "parent_id"),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(segment_length))



segment_counts <- spines_seg %>%
  group_by(mouse_id, cell_id, genotype, age, parent_id, depth, segment_length) %>%
  summarise(
    n_total = n(),
    n_thin = sum(spine_morphology == "Thin", na.rm = TRUE),
    n_stub = sum(spine_morphology == "Stubby", na.rm = TRUE),
    n_mush = sum(spine_morphology == "Mushroom", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    density_total = n_total / segment_length,
    density_thin = n_thin / segment_length,
    density_stub = n_stub / segment_length,
    density_mush = n_mush / segment_length
  )

all_results <- data.frame()

run_anova <- function(df, response, depth_value) {
  
  if (var(df[[response]], na.rm = TRUE) == 0) {
    return(data.frame(
      depth = depth_value,
      response = response,
      result_type = "skipped_no_variance"
    ))
  }
  
  f <- as.formula(paste(response, "~ genotype * age"))
  fit <- aov(f, data = df)
  
  # ANOVA results
  anova_df <- as.data.frame(summary(fit)[[1]])
  anova_df$Effect <- rownames(anova_df)
  
  anova_out <- anova_df %>%
    mutate(
      depth = depth_value,
      response = response,
      result_type = "ANOVA",
      contrast = NA,
      comparison_type = NA,
      p_value = `Pr(>F)`
    ) %>%
    select(depth, response, result_type, comparison_type, Effect, contrast, p_value, everything())
  
  em <- emmeans(fit, ~ genotype * age)
  
  # Sidak post hoc: WT vs KO within each age
  sidak_age <- as.data.frame(summary(pairs(em, by = "age", adjust = "sidak"))) %>%
    mutate(
      depth = depth_value,
      response = response,
      result_type = "Sidak_posthoc",
      comparison_type = "genotype_within_age",
      Effect = NA
    ) %>%
    rename(p_value = p.value) %>%
    select(depth, response, result_type, comparison_type, Effect, contrast, age, p_value, everything())
  
  # Sidak post hoc: 10w vs 20w within each genotype
  sidak_genotype <- as.data.frame(summary(pairs(em, by = "genotype", adjust = "sidak"))) %>%
    mutate(
      depth = depth_value,
      response = response,
      result_type = "Sidak_posthoc",
      comparison_type = "age_within_genotype",
      Effect = NA
    ) %>%
    rename(p_value = p.value) %>%
    select(depth, response, result_type, comparison_type, Effect, contrast, genotype, p_value, everything())
  
  bind_rows(anova_out, sidak_age, sidak_genotype)
}

depths <- sort(unique(segment_counts$depth))

for (d in depths) {
  df_d <- segment_counts %>% filter(depth == d)
  if (nrow(df_d) < 2) next
  
  all_results <- bind_rows(
    all_results,
    run_anova(df_d, "density_total", d),
    run_anova(df_d, "density_thin", d),
    run_anova(df_d, "density_mush", d),
    run_anova(df_d, "density_stub", d)
  )
}

significant_results <- all_results %>%
  filter(
    !is.na(p_value),
    p_value < 0.05
  )

output_dir <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/3_morphology/5_stat_summary"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_path <- file.path(
  output_dir,
  "spine_density_stat_summary.csv"
)

# Save ONLY the full results table
write_csv(all_results, output_path)

cat("Significant results only (p < 0.05):\n\n")

print(
  as.data.frame(significant_results),
  row.names = FALSE
)

cat("\nALL results saved to:\n", output_path, "\n")