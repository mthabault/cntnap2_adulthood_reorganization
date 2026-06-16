library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)


path <- #paste path here
df <- read.csv(path)

# Keep only 1-step transitions
df_1step <- df %>%
  filter(s1_len == 1, s2_len == 1)

behaviors <- c(
  "digging",
  "unsupported_rearing",
  "supported_rearing",
  "scratching",
  "head_body_twitch",
  "grooming_short",
  "grooming_long"
)


# Wide table
df_wide <- df_1step %>%
  mutate(
    s1 = factor(s1, levels = behaviors),
    s2 = factor(s2, levels = behaviors),
    pair = paste(s1, s2, sep = "→")
  ) %>%
  group_by(mouse, genotype, age, pair) %>%
  summarise(prob = mean(p_s2_given_s1), .groups = "drop") %>%
  pivot_wider(
    names_from = pair,
    values_from = prob,
    values_fill = 0
  )

pairs <- names(df_wide)[-(1:3)]


# Split groups
wt10 <- df_wide %>% filter(genotype == "wt", age == "10w")
ko10 <- df_wide %>% filter(genotype == "ko", age == "10w")
wt20 <- df_wide %>% filter(genotype == "wt", age == "20w")
ko20 <- df_wide %>% filter(genotype == "ko", age == "20w")

wt10_vals <- wt10 %>% select(all_of(pairs))
ko10_vals <- ko10 %>% select(all_of(pairs))
wt20_vals <- wt20 %>% select(all_of(pairs))
ko20_vals <- ko20 %>% select(all_of(pairs))


# Wilcoxon tests + FDR

# WT vs KO at 10w
res_wtko_10w <- lapply(pairs, function(p) {
  test <- wilcox.test(wt10_vals[[p]], ko10_vals[[p]], exact = FALSE)
  c(W = unname(test$statistic), p = test$p.value)
})
res_wtko_10w <- do.call(rbind, res_wtko_10w)
W_wtko_10w <- res_wtko_10w[, "W"]
pvals_wtko_10w <- res_wtko_10w[, "p"]
pvals_wtko_10w_fdr <- p.adjust(pvals_wtko_10w, method = "fdr")

# WT vs KO at 20w
res_wtko_20w <- lapply(pairs, function(p) {
  test <- wilcox.test(wt20_vals[[p]], ko20_vals[[p]], exact = FALSE)
  c(W = unname(test$statistic), p = test$p.value)
})
res_wtko_20w <- do.call(rbind, res_wtko_20w)
W_wtko_20w <- res_wtko_20w[, "W"]
pvals_wtko_20w <- res_wtko_20w[, "p"]
pvals_wtko_20w_fdr <- p.adjust(pvals_wtko_20w, method = "fdr")

# WT 10w vs WT 20
common_wt <- intersect(wt10$mouse, wt20$mouse)
wt10_paired <- wt10 %>% filter(mouse %in% common_wt) %>% arrange(mouse)
wt20_paired <- wt20 %>% filter(mouse %in% common_wt) %>% arrange(mouse)

res_wt_age <- lapply(pairs, function(p) {
  test <- wilcox.test(wt10_paired[[p]], wt20_paired[[p]], paired = TRUE, exact = FALSE)
  c(W = unname(test$statistic), p = test$p.value)
})
res_wt_age <- do.call(rbind, res_wt_age)
W_wt_age <- res_wt_age[, "W"]
pvals_wt_age <- res_wt_age[, "p"]
pvals_wt_age_fdr <- p.adjust(pvals_wt_age, method = "fdr")

# KO 10w vs KO 20w
common_ko <- intersect(ko10$mouse, ko20$mouse)
ko10_paired <- ko10 %>% filter(mouse %in% common_ko) %>% arrange(mouse)
ko20_paired <- ko20 %>% filter(mouse %in% common_ko) %>% arrange(mouse)

res_ko_age <- lapply(pairs, function(p) {
  test <- wilcox.test(ko10_paired[[p]], ko20_paired[[p]], paired = TRUE, exact = FALSE)
  c(W = unname(test$statistic), p = test$p.value)
})
res_ko_age <- do.call(rbind, res_ko_age)
W_ko_age <- res_ko_age[, "W"]
pvals_ko_age <- res_ko_age[, "p"]
pvals_ko_age_fdr <- p.adjust(pvals_ko_age, method = "fdr")

# Mean differences
diff_ko_wt_10w <- colMeans(ko10_vals) - colMeans(wt10_vals)
diff_ko_wt_20w <- colMeans(ko20_vals) - colMeans(wt20_vals)
diff_wt_age    <- colMeans(wt20_vals) - colMeans(wt10_vals)
diff_ko_age    <- colMeans(ko20_vals) - colMeans(ko10_vals)


# Long-format
diff_long <- function(diff_vec, contrast_name) {
  tibble(
    pair = pairs,
    diff = as.numeric(diff_vec)
  ) %>%
    separate(pair, into = c("s1", "s2"), sep = "→") %>%
    mutate(
      s1 = factor(s1, levels = behaviors),
      s2 = factor(s2, levels = behaviors),
      contrast = contrast_name
    )
}

p_long <- function(pvec_fdr, pvec_raw, wvec, contrast_name) {
  tibble(
    pair = pairs,
    p = as.numeric(pvec_raw),
    p_fdr = as.numeric(pvec_fdr),
    W = as.numeric(wvec)
  ) %>%
    separate(pair, into = c("s1", "s2"), sep = "→") %>%
    mutate(
      s1 = factor(s1, levels = behaviors),
      s2 = factor(s2, levels = behaviors),
      contrast = contrast_name
    )
}

# Plotting tables
df_mats <- bind_rows(
  diff_long(diff_ko_wt_10w, "KO − WT (10w)"),
  diff_long(diff_ko_wt_20w, "KO − WT (20w)"),
  diff_long(diff_wt_age,    "WT: 20w − 10w"),
  diff_long(diff_ko_age,    "KO: 20w − 10w")
) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c(
        "KO − WT (10w)",
        "KO − WT (20w)",
        "WT: 20w − 10w",
        "KO: 20w − 10w"
      )
    )
  )

df_pvals <- bind_rows(
  p_long(pvals_wtko_10w_fdr, pvals_wtko_10w, W_wtko_10w, "KO − WT (10w)"),
  p_long(pvals_wtko_20w_fdr, pvals_wtko_20w, W_wtko_20w, "KO − WT (20w)"),
  p_long(pvals_wt_age_fdr,   pvals_wt_age,   W_wt_age,   "WT: 20w − 10w"),
  p_long(pvals_ko_age_fdr,   pvals_ko_age,   W_ko_age,   "KO: 20w − 10w")
)

df_mats_sig <- df_mats %>%
  left_join(df_pvals, by = c("s1", "s2", "contrast")) %>%
  mutate(
    sig_label = case_when(
      is.na(p_fdr) ~ "",
      p_fdr < 1e-4 ~ "****",
      p_fdr < 1e-3 ~ "***",
      p_fdr < 1e-2 ~ "**",
      p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


# Print results table
full_table <- df_mats_sig %>%
  arrange(contrast, s1, s2) %>%
  select(contrast, s1, s2, diff, W, p, p_fdr, sig_label)

print(full_table, n = Inf)

# Plots

df_genotype_sig <- df_mats_sig %>%
  filter(contrast %in% c("KO − WT (10w)", "KO − WT (20w)"))

df_age_sig <- df_mats_sig %>%
  filter(contrast %in% c("WT: 20w − 10w", "KO: 20w − 10w"))

lim <- max(abs(df_mats$diff), na.rm = TRUE)


p_genotype <- ggplot(df_genotype_sig, aes(x = s1, y = s2, fill = diff)) +
  geom_tile(color = "white") +
  geom_text(
    data = filter(df_genotype_sig, sig_label != ""),
    aes(label = sig_label),
    size = 3,
    color = "black"
  ) +
  scale_fill_gradient2(
    name = "Δ P(s2 | s1)",
    midpoint = 0,
    limits = c(-lim, lim)
  ) +
  scale_y_discrete(limits = rev) +
  facet_wrap(~ contrast, nrow = 1) +
  labs(
    x = "State 1 (origin)",
    y = "State 2 (following)"
  ) +
  theme_prism(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(8, "pt")
  )

p_age <- ggplot(df_age_sig, aes(x = s1, y = s2, fill = diff)) +
  geom_tile(color = "white") +
  geom_text(
    data = filter(df_age_sig, sig_label != ""),
    aes(label = sig_label),
    size = 3,
    color = "black"
  ) +
  scale_fill_gradient2(
    name = "Δ P(s2 | s1)",
    midpoint = 0,
    limits = c(-lim, lim)
  ) +
  scale_y_discrete(limits = rev) +
  facet_wrap(~ contrast, nrow = 1) +
  labs(
    x = "State 1 (origin)",
    y = "State 2 (following)"
  ) +
  theme_prism(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(8, "pt")
  )

print(p_genotype)
print(p_age)

# Save stats

# Output path
out_path <- #paste path here

# Save statistics table
write.csv(full_table, out_path, row.names = FALSE)

cat("\nStat summary saved to:\n", out_path, "\n")