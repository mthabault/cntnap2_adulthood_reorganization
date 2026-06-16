library(dplyr)
library(tidyr)
library(readr)
library(vegan)

file_path <- #paste path here
df <- read_csv(file_path)

# Prepare data
# Define behaviours
exploratory <- c("digging", "supported_rearing", "unsupported_rearing")
stereotyped <- c("scratching", "head_body_twitch", "grooming_long", "grooming_short")
all_beh <- c(exploratory, stereotyped)
df_T <- df %>%    #store the total duration of each video for further analysis
  group_by(video) %>%
  summarise(T = max(end_time), .groups = "drop")

# Compute activity and bias values for behavioural space
df_sub <- df %>%
  filter(behavior %in% all_beh) %>%
  group_by(video, behavior) %>%
  summarise(dur = sum(duration), .groups = "drop") %>%
  pivot_wider(
    names_from = behavior,
    values_from = dur,
    values_fill = 0
  )

df_unit <- df_sub %>%
  left_join(df_T, by = "video")

df_unit <- df_unit %>%
  mutate(
    dur_explor = digging + supported_rearing + unsupported_rearing,
    dur_stereo = scratching + head_body_twitch + grooming_long + grooming_short,
    dur_total  = dur_explor + dur_stereo
  )

df_unit <- df_unit %>%
  mutate(
    Bias = if_else(
      dur_total > 0,
      (dur_explor - dur_stereo) / dur_total,   # -1 (pure stereotyped) → +1 (pure exploratory)
      0
    ),
    Activity = dur_total / T                   # fraction of session with scored behaviour
  )


# Add metadata in dataset
df_meta <- df %>%
  select(video, mouse, genotype, sex, age) %>%
  distinct()
df_unit <- df_unit %>%
  left_join(df_meta, by = "video") %>%
  mutate(group = paste(genotype, age, sep = "_"))

# Create vectors
# Sanity check = make sure all animals have data at both ages
paired <- df_unit %>%
  group_by(mouse) %>%
  filter(n_distinct(age) == 2) %>%
  ungroup()

arrows_mouse <- paired %>%
  select(mouse, genotype, age, Bias, Activity) %>%
  arrange(mouse, age) %>%   # makes sure 10w comes before 20w
  group_by(mouse, genotype) %>%
  summarise(
    x    = Bias[age == "10w"],
    y    = Activity[age == "10w"],
    xend = Bias[age == "20w"],
    yend = Activity[age == "20w"],
    .groups = "drop"
  )
print(arrows_mouse)
print(nrow(arrows_mouse))  # This is the end of the sanity check = check how many pairs have been included.
# should match the sample size.


#Analysis

# Set vectors 
vecs <- arrows_mouse %>%
  mutate(
    genotype = as.factor(trimws(tolower(genotype))),  # ensure "ko"/"wt"
    dx       = xend - x,
    dy       = yend - y
  )

print(vecs)
print(table(vecs$genotype))

# dx/dy correlations 
cor_summary <- vecs %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    cor_test = list(cor.test(dx, dy, method = "pearson")),
    cor_dx_dy = cor_test[[1]]$estimate,
    ci_low    = cor_test[[1]]$conf.int[1],
    ci_high   = cor_test[[1]]$conf.int[2],
    p_value   = cor_test[[1]]$p.value,
    .groups = "drop"
  )
print(cor_summary)

# Fisher r-to-z test for difference in correlations
g1_n  <- cor_summary$n[1]
g2_n  <- cor_summary$n[2]
g1_r  <- cor_summary$cor_dx_dy[1]
g2_r  <- cor_summary$cor_dx_dy[2]

z1 <- 0.5 * log((1 + g1_r) / (1 - g1_r))
z2 <- 0.5 * log((1 + g2_r) / (1 - g2_r))

z_diff <- (z1 - z2) / sqrt(1 / (g1_n - 3) + 1 / (g2_n - 3))
p_diff <- 2 * pnorm(-abs(z_diff))

cat("\nFisher r-to-z test for difference in dx–dy correlation between genotypes:\n")
cat("z =", z_diff, "  p =", p_diff, "\n")


# Procrustes / Protest
ko_mat <- as.matrix(vecs %>% filter(genotype == gen_levels[1]) %>% select(dx, dy))
wt_mat <- as.matrix(vecs %>% filter(genotype == gen_levels[2]) %>% select(dx, dy))

if (nrow(ko_mat) == nrow(wt_mat)) {
  proc_res   <- procrustes(wt_mat, ko_mat, scale = TRUE)
  protest_res <- protest(wt_mat, ko_mat, permutations = 999)
  
  cat("\nProcrustes analysis (configuration similarity KO vs WT in shift space):\n")
  print(proc_res)
  cat("\nProtest permutation test:\n")
  print(protest_res)
} else {
  cat("\nProcrustes/Protest not run: unequal sample sizes between genotypes.\n")
}


# Save all stats in csv file

# Table for correlations
cor_out <- cor_summary %>%
  select(genotype, n, cor_dx_dy, ci_low, ci_high, p_value)

# Table for Fisher test
fisher_out <- tibble(
  test = "fisher_r_to_z",
  genotype_1 = cor_summary$genotype[1],
  genotype_2 = cor_summary$genotype[2],
  n1 = g1_n,
  n2 = g2_n,
  r1 = g1_r,
  r2 = g2_r,
  z1 = z1,
  z2 = z2,
  z_diff = z_diff,
  p_value = p_diff
)

# Procrustes / Protest results
if (exists("proc_res") & exists("protest_res")) {
  
  proc_out <- tibble(
    test = "procrustes",
    ss = proc_res$ss
  )
  
  protest_out <- tibble(
    test = "protest",
    t0 = protest_res$t0,
    p_value = protest_res$signif
  )
  
} else {
  proc_out <- tibble()
  protest_out <- tibble()
}

# Combine everything into one dataframe
final_stats <- bind_rows(
  cor_out %>% mutate(test = "correlation"),
  fisher_out,
  proc_out,
  protest_out
)

# Define output path
out_path <- #past path here

# Save file
write_csv(final_stats, out_path)

cat("\nStat summary saved to:\n", out_path, "\n")