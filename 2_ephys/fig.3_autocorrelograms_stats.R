library(tidyverse)
library(emmeans)

# Load data
data <- read.csv(#paste path here)

# Prepare variables
data <- data %>%
  mutate(
    genotype  = factor(genotype),
    age       = factor(age),
    group     = interaction(genotype, age, sep = "_"),
    cell_id   = factor(recording_id),
    lag_ms    = as.numeric(lag_ms)
  )

# Compute bin width
bin_width_ms <- data %>%
  distinct(lag_ms) %>%
  arrange(lag_ms) %>%
  mutate(diff = lag(lag_ms)) %>%
  filter(!is.na(diff), diff > 0) %>%
  summarise(bw = median(diff)) %>%
  pull(bw)

cat("Bin width (ms):", bin_width_ms, "\n\n")

# Compute autocorrelation shape descriptors for 0–500 ms
descriptors <- data %>%
  filter(lag_ms >= 0, lag_ms <= 500) %>%
  group_by(cell_id, genotype, age, group) %>%
  arrange(lag_ms) %>%
  summarise(
    auc_0_100 = {
      x <- lag_ms[lag_ms <= 100]
      y <- autocorr_value[lag_ms <= 100]
      if (length(x) < 2) NA_real_ else sum(y) * bin_width_ms
    },
    auc_100_500 = {
      x <- lag_ms[lag_ms > 100 & lag_ms <= 500]
      y <- autocorr_value[lag_ms > 100 & lag_ms <= 500]
      if (length(x) < 2) NA_real_ else sum(y) * bin_width_ms
    },
    peak_width_ms = {
      x <- lag_ms[lag_ms <= 200]
      y <- autocorr_value[lag_ms <= 200]
      if (length(x) < 3) {
        NA_real_
      } else {
        max_y <- max(y, na.rm = TRUE)
        half  <- max_y / 2
        idx   <- which(y >= half)
        if (length(idx) == 0) {
          NA_real_
        } else {
          x[max(idx)] - x[min(idx)]
        }
      }
    },
    decay_slope = {
      x <- lag_ms[lag_ms >= 50 & lag_ms <= 500]
      y <- autocorr_value[lag_ms >= 50 & lag_ms <= 500]
      if (length(x) < 3) {
        NA_real_
      } else {
        coef(lm(y ~ x))[2]
      }
    },
    sd_autocorr = sd(autocorr_value, na.rm = TRUE),
    .groups = "drop"
  )



cat("Head descriptors:\n")
print(head(descriptors))
cat("\n")

# Helper function to fit LM and emmeans for one descriptor
analyze_descriptor <- function(df, response) {
  formula_ga <- as.formula(paste(response, "~ genotype * age"))
  formula_grp <- as.formula(paste(response, "~ group"))
  
  cat("==================================================\n")
  cat("Descriptor:", response, "\n\n")
  
  cat("Linear model: ", deparse(formula_ga), "\n")
  model_ga <- lm(formula_ga, data = df)
  print(anova(model_ga))
  cat("\n")
  
  cat("Linear model with group (for pairwise comparisons): ", deparse(formula_grp), "\n")
  model_grp <- lm(formula_grp, data = df)
  emm <- emmeans(model_grp, pairwise ~ group, adjust = "tukey")
  
  cat("\nEstimated marginal means by group:\n")
  print(emm$emmeans)
  cat("\nPairwise Tukey comparisons:\n")
  print(emm$contrasts)
  cat("\n")
}

# Run analyses for each descriptor
analyze_descriptor(descriptors, "auc_0_100")
analyze_descriptor(descriptors, "auc_100_500")
analyze_descriptor(descriptors, "peak_width_ms")
analyze_descriptor(descriptors, "decay_slope")
analyze_descriptor(descriptors, "sd_autocorr")
