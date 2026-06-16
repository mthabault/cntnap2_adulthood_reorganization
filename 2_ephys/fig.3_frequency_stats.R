# Load required libraries
library(dplyr)

# Load the data
data <- read.csv(#paste path here)

data$group <- interaction(data$genotype, data$age, drop = TRUE)

# t-tests
t_results <- lapply(c("10w", "20w"), function(age_group) {
  df <- filter(data, age == age_group)
  t_out <- t.test(frequency ~ genotype, data = df, var.equal = TRUE)
  
  data.frame(
    test = "t_test_WT_vs_KO",
    age = age_group,
    t_statistic = unname(t_out$statistic),
    df = unname(t_out$parameter),
    p_value = t_out$p.value,
    mean_wt = mean(df$frequency[df$genotype == "wt"]),
    mean_ko = mean(df$frequency[df$genotype == "ko"])
  )
}) %>% bind_rows()

# Save t-test results
write.csv(
  t_results,
  "#path here",
  row.names = FALSE
)

cat("All statistical results saved successfully.\n")