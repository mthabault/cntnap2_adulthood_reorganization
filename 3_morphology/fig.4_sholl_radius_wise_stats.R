library(tidyverse)
library(broom)

sholl_data <- read.csv(
  "#path"
)

sholl_data <- sholl_data %>%
  mutate(
    mouse_id      = factor(mouse_id),
    cell_id       = factor(cell_id),
    genotype      = factor(genotype),
    age           = factor(age),
    radius        = as.numeric(radius),
    intersections = as.numeric(intersections),
    cell_uid      = interaction(mouse_id, cell_id, drop = TRUE)
  )

sholl_data$genotype <- relevel(sholl_data$genotype, ref = "wt")

analyze_age_per_radius <- function(dat, age_label) {
  message("")
  message(paste("per-radius analysis for age =", age_label, "(ko vs wt)"))
  message("")
  
  dat_age <- dat %>%
    filter(age == age_label, genotype %in% c("wt", "ko")) %>%
    droplevels()
  
  lm_results <- dat_age %>%
    group_by(radius) %>%
    group_modify(~ {
      fit <- lm(intersections ~ genotype, data = .x)
      coef_tbl <- tidy(fit)
      info_tbl <- glance(fit)
      term_row <- coef_tbl %>% filter(term == "genotypeko")
      tibble(
        radius      = unique(.x$radius),
        estimate    = term_row$estimate,
        std.error   = term_row$std.error,
        t.value     = term_row$statistic,
        p.value_raw = term_row$p.value,
        df          = info_tbl$df.residual
      )
    }) %>%
    ungroup()
  
  k <- nrow(lm_results)
  lm_results <- lm_results %>%
    mutate(
      p.fdr    = p.adjust(p.value_raw, method = "fdr"),
      age      = age_label,
      contrast = "ko_vs_wt"
    )
  
  print(lm_results, n = Inf)
  
  curve_summary <- dat_age %>%
    group_by(genotype, radius) %>%
    summarise(
      mean_intersections = mean(intersections, na.rm = TRUE),
      se_intersections   = sd(intersections, na.rm = TRUE) / sqrt(n()),
      n                  = n(),
      .groups = "drop"
    )
  
  print(curve_summary)
  
  invisible(list(
    data_age       = dat_age,
    lm_results     = lm_results,
    curve_summary  = curve_summary
  ))
}

# Per-radius comparison of ages within a genotype
analyze_genotype_per_radius <- function(dat, genotype_label) {
  message("")
  message(paste("per-radius analysis for genotype =", genotype_label, "(20w vs 10w)"))
  message("")
  
  dat_gen <- dat %>%
    filter(genotype == genotype_label, age %in% c("10w", "20w")) %>%
    droplevels()
  
  # ensure age is a factor with 10w as reference
  dat_gen <- dat_gen %>%
    mutate(age = factor(age))
  dat_gen$age <- relevel(dat_gen$age, ref = "10w")
  
  lm_results <- dat_gen %>%
    group_by(radius) %>%
    group_modify(~ {
      fit <- lm(intersections ~ age, data = .x)
      coef_tbl <- tidy(fit)
      info_tbl <- glance(fit)
      term_row <- coef_tbl %>% filter(term == "age20w")
      tibble(
        radius      = unique(.x$radius),
        estimate    = term_row$estimate,
        std.error   = term_row$std.error,
        t.value     = term_row$statistic,
        p.value_raw = term_row$p.value,
        df          = info_tbl$df.residual
      )
    }) %>%
    ungroup()
  
  k <- nrow(lm_results)
  lm_results <- lm_results %>%
    mutate(
      p.fdr    = p.adjust(p.value_raw, method = "fdr"),
      age      = "20w_vs_10w",
      contrast = paste0("age_20w_vs_10w_in_", genotype_label)
    )
  
  print(lm_results)
  
  curve_summary <- dat_gen %>%
    group_by(age, radius) %>%
    summarise(
      mean_intersections = mean(intersections, na.rm = TRUE),
      se_intersections   = sd(intersections, na.rm = TRUE) / sqrt(n()),
      n                  = n(),
      .groups = "drop"
    )
  
  print(curve_summary)
  
  invisible(list(
    data_gen       = dat_gen,
    lm_results     = lm_results,
    curve_summary  = curve_summary
  ))
}

res_10w <- analyze_age_per_radius(sholl_data, "10w")   # ko vs wt at 10w
res_20w <- analyze_age_per_radius(sholl_data, "20w")   # ko vs wt at 20w

res_wt <- analyze_genotype_per_radius(sholl_data, "wt")  # 20w vs 10w in wt
res_ko <- analyze_genotype_per_radius(sholl_data, "ko")  # 20w vs 10w in ko

results_fdr <- bind_rows(
  res_10w$lm_results,
  res_20w$lm_results,
  res_wt$lm_results,
  res_ko$lm_results
)

write.csv(
  results_fdr,
  "#path",
  row.names = FALSE
)
