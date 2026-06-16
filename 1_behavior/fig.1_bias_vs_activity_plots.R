library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)
library(ggforce)
library(readr)
library (Cairo)

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


# Create plots
  # Define colors f
colors <- c(
  "wt_10w" = "grey",
  "ko_10w" = "#34d6fa",
  "wt_20w" = "black",
  "ko_20w" = "#0213f7"
)

  # Plot function
behaviour_plot <- function(data, title) {
  ggplot(data, aes(x = Bias, y = Activity, color = group)) +
    geom_point(size = 3, alpha = 1) +
    stat_ellipse(type = "norm", level = 0.95, linewidth = 1.2) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_x_continuous(limits = c(-1.5, 1.5), breaks = seq(-1, 1, 0.5)) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_prism() +
    labs(
      title = title,
      x = "Bias (exploratory →   /   stereotyped ←)",
      y = "Activity (fraction of session with scored behaviour)",
      color = "Group"
    )
}

# WT Plot
df_wt <- df_unit %>% filter(genotype == "wt")
p_wt_age <- behaviour_plot(df_wt, "WT: 10w vs 20w") +
  geom_segment(
    data = arrows_mouse %>% filter(genotype == "wt"),
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.02, "npc")),
    inherit.aes = FALSE,
    linewidth = 0.4,
    alpha = 1,
    color = "black"
  )

# KO plot
df_ko <- df_unit %>% filter(genotype == "ko")
p_ko_age <- behaviour_plot(df_ko, "KO: 10w vs 20w") +
  geom_segment(
    data = arrows_mouse %>% filter(genotype == "ko"),
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.02, "npc")),
    inherit.aes = FALSE,
    linewidth = 0.4,
    alpha = 1,
    color = "black"
  )

print(p_wt_age)
print(p_ko_age)


# 12. Save the plots as EPS
output_dir <- "paste your path here"

plots <- list(
  WT_10w_vs_20w = p_wt_age,
  KO_10w_vs_20w = p_ko_age
)

for (nm in names(plots)) {
  cat("\n\n===== Saving ", nm, " =====\n")
  p <- plots[[nm]]
  print(p)
  
  eps_filename <- file.path(output_dir, paste0(nm, ".eps"))
  
  ggsave(
    filename = eps_filename,
    plot = p,
    device = cairo_ps,
    width = 8,
    height = 8,
    dpi = 300,
    bg = "white",
    fallback_resolution = 300,
    antialias = "none"
  )
}