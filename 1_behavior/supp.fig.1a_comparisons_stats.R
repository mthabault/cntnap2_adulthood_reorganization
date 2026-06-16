library(dplyr)
library(readr)
library(afex)
library(emmeans)

# make afex play nice with emmeans / type III
afex_options(type = 3)
options(contrasts = c("contr.sum", "contr.poly"))

behaviour_data <- read.csv(
  "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/1_behaviour/4_csv_datasets/behaviour_summary_per_mouse.csv",
  stringsAsFactors = FALSE
)

composites <- behaviour_data %>%
  mutate(
    composite = case_when(
      behavior %in% c(
        "scratching",
        "head_body_twitch",
        "grooming_long",
        "grooming_short"
      ) ~ "repetitive and stereotyped behaviour",
      
      behavior %in% c(
        "digging",
        "supported_rearing",
        "unsupported_rearing"
      ) ~ "exploratory behaviour",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(composite)) %>%
  group_by(video, mouse, genotype, age, composite) %>%
  summarise(
    n_behaviours = sum(n_behaviours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(behavior = composite)

behaviour_data2 <- bind_rows(
  behaviour_data,
  composites
)

behaviours <- c(
  "repetitive and stereotyped behaviour",
  "exploratory behaviour"
)

behaviour_data2 <- behaviour_data2 %>%
  filter(
    behavior %in% behaviours,
    age %in% c("10w", "20w")
  ) %>%
  mutate(
    genotype = tolower(genotype),
    genotype = factor(genotype, levels = c("wt", "ko")),
    age      = factor(age, levels = c("10w", "20w")),
    mouse    = factor(mouse)
  )

all_results <- data.frame()

run_rmanova_counts <- function(df) {
  
  behaviour_name <- unique(df$behavior)
  
  cat("\nBehaviour:", behaviour_name, "\n\n")
  
  # rmANOVA
  fit <- aov_ez(
    id = "mouse",
    dv = "n_behaviours",
    data = df,
    within = "age",
    between = "genotype",
    type = 3,
    factorize = FALSE
  )
  
  # ANOVA table
  anova_df <- as.data.frame(fit$anova_table)
  anova_df$Effect <- rownames(anova_df)
  
  anova_out <- anova_df %>%
    mutate(
      behavior = behaviour_name,
      result_type = "RM_ANOVA"
    ) %>%
    rename(
      p_value = `Pr(>F)`
    ) %>%
    select(
      behavior,
      result_type,
      Effect,
      everything()
    )
  
  # Significant ANOVA only
  sig_anova <- anova_out %>%
    filter(
      !is.na(p_value),
      p_value < 0.05
    )
  
  cat("Significant RM-ANOVA effects:\n")
  print(sig_anova)
  
  # genotype differences at each age
  emm_geno_by_age <- emmeans(
    fit,
    ~ genotype | age
  )
  
  geno_posthoc <- as.data.frame(
    summary(
      pairs(
        emm_geno_by_age,
        adjust = "tukey"
      )
    )
  ) %>%
    mutate(
      behavior = behaviour_name,
      result_type = "Tukey_posthoc",
      comparison_type = "genotype_within_age"
    ) %>%
    rename(
      p_value = p.value
    )
  
  # age differences within genotype
  emm_age_by_geno <- emmeans(
    fit,
    ~ age | genotype
  )
  
  age_posthoc <- as.data.frame(
    summary(
      pairs(
        emm_age_by_geno,
        adjust = "tukey"
      )
    )
  ) %>%
    mutate(
      behavior = behaviour_name,
      result_type = "Tukey_posthoc",
      comparison_type = "age_within_genotype"
    ) %>%
    rename(
      p_value = p.value
    )
  
  # Significant posthoc only
  sig_posthoc <- bind_rows(
    geno_posthoc,
    age_posthoc
  ) %>%
    filter(
      !is.na(p_value),
      p_value < 0.05
    )
  
  cat("\nSignificant Tukey posthoc results:\n")
  print(sig_posthoc)
  
  # Return all results
  bind_rows(
    anova_out,
    geno_posthoc,
    age_posthoc
  )
}

for (beh in unique(behaviour_data2$behavior)) {
  
  df_beh <- behaviour_data2 %>%
    filter(behavior == beh)
  
  results_beh <- run_rmanova_counts(df_beh)
  
  all_results <- bind_rows(
    all_results,
    results_beh
  )
}

# Output directory
output_dir <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/1_behaviour/6_stat_summary"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Output file
output_path <- file.path(
  output_dir,
  "supp_fig1a_stat_summary.csv"
)

# Save ALL results
write_csv(
  all_results,
  output_path
)

cat("\nALL results saved to:\n")
cat(output_path, "\n")