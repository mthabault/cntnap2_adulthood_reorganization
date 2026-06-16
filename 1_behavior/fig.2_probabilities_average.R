library(dplyr)
library(tidyr)
library(emmeans)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggprism)

input_path <- #past path here
df <- read.csv(input_path)

behaviours <- c(
  "digging",
  "supported_rearing",
  "unsupported_rearing",
  "scratching",
  "grooming_long",
  "grooming_short",
  "head_body_twitch"
)

empty <- "empty"

# Transition classifier
classify_transition <- function(b1_row, next_rows) {
  t1_end <- b1_row$end_time
  b1     <- b1_row$behavior
  
  if (nrow(next_rows) == 0) return(NULL)
  
  for (i in seq_len(nrow(next_rows))) {
    beh2   <- next_rows$behavior[i]
    start2 <- next_rows$start_time[i]
    delta  <- start2 - t1_end
    
    if (beh2 == empty) {
      if (next_rows$duration[i] >= 3) return("isolated")
      next
    }
    
    if (delta >= 3) return("isolated")
    
    if (beh2 == b1) return("same") else return("diff")
  }
  
  "isolated"
}

# Order events and compute transitions
df <- df %>%
  arrange(mouse, video, start_time) %>%
  group_by(mouse, genotype, sex, age, video) %>%
  mutate(idx = dplyr::row_number()) %>%
  ungroup()

groups <- df %>%
  group_by(mouse, genotype, sex, age, video) %>%
  group_split()

trans_list <- list()

for (g in seq_along(groups)) {
  subdf <- groups[[g]]
  origin <- subdf[subdf$behavior %in% behaviours, ]
  if (nrow(origin) < 2) next
  
  for (i in 1:(nrow(origin) - 1)) {
    row1 <- origin[i, ]
    next_rows <- subdf[subdf$idx > row1$idx, ]
    cat_type <- classify_transition(row1, next_rows)
    if (is.null(cat_type)) next
    
    trans_list[[length(trans_list) + 1]] <- data.frame(
      mouse    = row1$mouse,
      genotype = row1$genotype,
      sex      = row1$sex,
      age      = row1$age,
      behavior = row1$behavior,
      trans    = cat_type,
      stringsAsFactors = FALSE
    )
  }
}

df_trans <- dplyr::bind_rows(trans_list)

# Counts and probabilities per mouse x age (averaged over behaviours)
df_counts <- df_trans %>%
  group_by(mouse, genotype, sex, age, behavior, trans) %>%
  tally(name = "n") %>%
  ungroup() %>%
  tidyr::pivot_wider(
    names_from  = trans,
    values_from = n,
    values_fill = 0
  ) %>%
  rename(
    n_same = same,
    n_diff = diff,
    n_iso  = isolated
  ) %>%
  mutate(
    n_events   = n_same + n_diff + n_iso,
    p_same     = n_same / n_events,
    p_diff     = n_diff / n_events,
    p_isolated = n_iso  / n_events
  )

df_mouse_age <- df_counts %>%
  group_by(mouse, genotype, age) %>%
  summarise(
    p_same     = mean(p_same, na.rm = TRUE),
    p_diff     = mean(p_diff, na.rm = TRUE),
    p_isolated = mean(p_isolated, na.rm = TRUE),
    .groups = "drop"
  )

# Sanity check: keep only mice with both ages
df_pairs <- df_mouse_age %>%
  filter(age %in% c("10w", "20w")) %>%
  group_by(mouse) %>%
  filter(n_distinct(age) == 2) %>%
  ungroup() %>%
  mutate(
    genotype = factor(genotype, levels = c("wt", "ko")),
    age      = factor(age,      levels = c("10w", "20w"))
  )

df_wide <- df_pairs %>%
  tidyr::pivot_wider(
    id_cols = c(mouse, genotype),
    names_from = age,
    values_from = c(p_same, p_diff, p_isolated)
  )

# Build long data for lmm & run llm
df_long_lmm <- df_wide %>%
  pivot_longer(
    cols = c(p_same_10w, p_same_20w,
             p_diff_10w, p_diff_20w,
             p_isolated_10w, p_isolated_20w),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(
    metric = case_when(
      grepl("p_same", var)     ~ "p_same",
      grepl("p_diff", var)     ~ "p_diff",
      grepl("p_isolated", var) ~ "p_isolated"
    ),
    age = ifelse(grepl("10w", var), "10w", "20w"),
    age = factor(age, levels = c("10w", "20w")),
    genotype = factor(genotype, levels = c("wt", "ko"))
  )

fit_metric <- function(df_long_lmm, metric_name) {
  sub <- df_long_lmm %>% filter(metric == metric_name)
  m <- lmer(value ~ genotype * age + (1 | mouse), data = sub)
  a <- anova(m, type = 3)
  em <- emmeans(m, pairwise ~ genotype * age, adjust = "tukey")
  list(
    data = sub,
    model = m,
    anova = a,
    emmeans = em
  )
}

cat("lmm for p_same\n")
res_same <- fit_metric(df_long_lmm, "p_same")
print(res_same$anova)
print(res_same$emmeans)

cat("\nlmm for p_diff\n")
res_diff <- fit_metric(df_long_lmm, "p_diff")
print(res_diff$anova)
print(res_diff$emmeans)

cat("\nlmm for p_isolated\n")
res_iso <- fit_metric(df_long_lmm, "p_isolated")
print(res_iso$anova)
print(res_iso$emmeans)

extract_p <- function(aov_tab, effect) {
  aov_tab[effect, "Pr(>F)"]
}

p_same_geno <- extract_p(res_same$anova, "genotype")
p_same_age  <- extract_p(res_same$anova, "age")
p_same_int  <- extract_p(res_same$anova, "genotype:age")

p_diff_geno <- extract_p(res_diff$anova, "genotype")
p_diff_age  <- extract_p(res_diff$anova, "age")
p_diff_int  <- extract_p(res_diff$anova, "genotype:age")

p_iso_geno <- extract_p(res_iso$anova, "genotype")
p_iso_age  <- extract_p(res_iso$anova, "age")
p_iso_int  <- extract_p(res_iso$anova, "genotype:age")

sig_report <- tibble(
  metric = c("p_same", "p_diff", "p_isolated"),
  p_genotype = c(p_same_geno, p_diff_geno, p_iso_geno),
  p_age = c(p_same_age, p_diff_age, p_iso_age),
  p_interaction = c(p_same_int, p_diff_int, p_iso_int)
)

print(sig_report)

# Plots
df_stack <- df_pairs %>%
  pivot_longer(
    cols = c(p_same, p_diff, p_isolated),
    names_to = "metric",
    values_to = "prob"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("p_same", "p_diff", "p_isolated"),
      labels = c("same", "diff", "isolated")
    ),
    mouse = factor(mouse)
  )

metric_colors <- c(
  same     = "#1b9e77",
  diff     = "#7570b3",
  isolated = "#d95f02"
)

# WT 
df_wt <- df_stack %>% filter(genotype == "wt")

gg_stack_wt <- ggplot(
  df_wt,
  aes(x = mouse, y = prob * 100, fill = metric)
) +
  geom_col() +
  facet_wrap(~ age, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = metric_colors) +
  labs(
    title = "WT",
    x = "mouse",
    y = "transition probability [%]",
    fill = "transition type"
  ) +
  theme_prism(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing.x = unit(10, "pt")
  )

print(gg_stack_wt)

# KO 
df_ko <- df_stack %>% filter(genotype == "ko")

gg_stack_ko <- ggplot(
  df_ko,
  aes(x = mouse, y = prob * 100, fill = metric)
) +
  geom_col() +
  facet_wrap(~ age, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = metric_colors) +
  labs(
    title = "KO",
    x = "mouse",
    y = "transition probability [%]",
    fill = "transition type"
  ) +
  theme_prism(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing.x = unit(10, "pt")
  )

print(gg_stack_ko)

# Compute mean values for each group
df_group_means <- df_stack %>%
  group_by(genotype, age, metric) %>%
  summarise(mean_prob = mean(prob), .groups = "drop") %>%
  mutate(mean_prob = mean_prob * 100)

# Helper function to make a stacked bar for a single group
make_mean_plot <- function(data, title_text) {
  ggplot(data, aes(x = "", y = mean_prob, fill = metric)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = metric_colors) +
    labs(
      title = title_text,
      x = NULL,
      y = "mean transition probability [%]",
      fill = "transition type"
    ) +
    theme_prism(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# WT 10w
plot_wt_10w <- make_mean_plot(
  df_group_means %>% filter(genotype == "wt", age == "10w"),
  "WT – 10w"
)

# WT 20w
plot_wt_20w <- make_mean_plot(
  df_group_means %>% filter(genotype == "wt", age == "20w"),
  "WT – 20w"
)

# KO 10w
plot_ko_10w <- make_mean_plot(
  df_group_means %>% filter(genotype == "ko", age == "10w"),
  "KO – 10w"
)

# KO 20w
plot_ko_20w <- make_mean_plot(
  df_group_means %>% filter(genotype == "ko", age == "20w"),
  "KO – 20w"
)


# Save all plots as EPS

output_dir <- "define your output directory"

# Mean plots (3 × 8)
mean_plots <- list(
  WT_10w = plot_wt_10w,
  WT_20w = plot_wt_20w,
  KO_10w = plot_ko_10w,
  KO_20w = plot_ko_20w
)

for (nm in names(mean_plots)) {
  
  eps_filename <- file.path(output_dir, paste0(nm, ".eps"))
  
  ggsave(
    filename = eps_filename,
    plot = mean_plots[[nm]],
    device = cairo_ps,
    width = 3,
    height = 8,
    dpi = 300
  )
}

# Stacked mouse plots (10 × 8)
stack_plots <- list(
  WT_stack = gg_stack_wt,
  KO_stack = gg_stack_ko
)

for (nm in names(stack_plots)) {
  
  eps_filename <- file.path(output_dir, paste0(nm, ".eps"))
  
  ggsave(
    filename = eps_filename,
    plot = stack_plots[[nm]],
    device = cairo_ps,
    width = 10,
    height = 8,
    dpi = 300
  )
}


# Save stats

# Output path
out_path <- "past path here"

# ANOVA tables
anova_same_out <- as.data.frame(res_same$anova) %>%
  tibble::rownames_to_column("effect") %>%
  mutate(metric = "p_same", test = "lmm_type3_anova")

anova_diff_out <- as.data.frame(res_diff$anova) %>%
  tibble::rownames_to_column("effect") %>%
  mutate(metric = "p_diff", test = "lmm_type3_anova")

anova_iso_out <- as.data.frame(res_iso$anova) %>%
  tibble::rownames_to_column("effect") %>%
  mutate(metric = "p_isolated", test = "lmm_type3_anova")

anova_out <- bind_rows(
  anova_same_out,
  anova_diff_out,
  anova_iso_out
)

# Estimated marginal means
emm_same_out <- as.data.frame(res_same$emmeans$emmeans) %>%
  mutate(metric = "p_same", test = "emmeans")

emm_diff_out <- as.data.frame(res_diff$emmeans$emmeans) %>%
  mutate(metric = "p_diff", test = "emmeans")

emm_iso_out <- as.data.frame(res_iso$emmeans$emmeans) %>%
  mutate(metric = "p_isolated", test = "emmeans")

emm_out <- bind_rows(
  emm_same_out,
  emm_diff_out,
  emm_iso_out
)

# Pairwise contrasts
contrast_same_out <- as.data.frame(res_same$emmeans$contrasts) %>%
  mutate(metric = "p_same", test = "tukey_pairwise_contrasts")

contrast_diff_out <- as.data.frame(res_diff$emmeans$contrasts) %>%
  mutate(metric = "p_diff", test = "tukey_pairwise_contrasts")

contrast_iso_out <- as.data.frame(res_iso$emmeans$contrasts) %>%
  mutate(metric = "p_isolated", test = "tukey_pairwise_contrasts")

contrast_out <- bind_rows(
  contrast_same_out,
  contrast_diff_out,
  contrast_iso_out
)

sig_report_out <- sig_report %>%
  mutate(test = "p_value_summary")

group_means_out <- df_group_means %>%
  mutate(test = "group_means_percent")

# Combine everything
final_stats <- bind_rows(
  anova_out,
  emm_out,
  contrast_out,
  sig_report_out,
  group_means_out
)

# Save CSV
write.csv(final_stats, out_path, row.names = FALSE)

cat("\nStat summary saved to:\n", out_path, "\n")