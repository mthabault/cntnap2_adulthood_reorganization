library(dplyr)
library(readr)

file_metrics <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/3_morphology/3_csv_datasets/sholl_metrics_per_cell.csv"

metrics <- read_csv(file_metrics) %>%
  mutate(
    genotype = factor(genotype),
    age      = factor(age),
    group    = factor(paste(genotype, age, sep = "_"))
  )

metric_names <- c(
  "AUC",
  "Imax",
  "Rmax",
  "field_radius",
  "critical_radius",
  "slope",
  "intercept",
  "prox_dist_ratio"
)

p_to_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")
           )
    )
  )
}

ages <- c("10w", "20w")

all_results <- data.frame()

for (m in metric_names) {
  
  cat("\n\nMetric:", m, "\n")
  
  for (a in ages) {
    
    cat("\nAge:", a, "\n")
    
    sub <- metrics %>%
      filter(age == a)
    
    # Welch t-test
    t_fit <- t.test(
      sub[[m]] ~ sub$genotype,
      var.equal = FALSE
    )
    
    p_t <- t_fit$p.value
    star_t <- p_to_stars(p_t)
    
    cat(
      "Welch t-test p =",
      p_t,
      star_t,
      "\n"
    )
    
    # Store results
    result_row <- data.frame(
      metric = m,
      age = a,
      test = "Welch_t_test",
      statistic_t = unname(t_fit$statistic),
      df = unname(t_fit$parameter),
      p_value = p_t,
      significance = star_t,
      mean_ko = mean(sub[[m]][sub$genotype == "ko"], na.rm = TRUE),
      mean_wt = mean(sub[[m]][sub$genotype == "wt"], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    all_results <- bind_rows(
      all_results,
      result_row
    )
  }
}

# Output directory
output_dir <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/3_morphology/5_stat_summary"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Output file
output_path <- file.path(
  output_dir,
  "sholl_metrics_stat_summary.csv"
)

# Save all results
write_csv(
  all_results,
  output_path
)

cat("\nALL results saved to:\n")
cat(output_path, "\n")

cat("\nDone.\n")