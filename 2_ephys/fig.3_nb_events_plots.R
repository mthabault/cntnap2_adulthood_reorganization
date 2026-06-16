library(tidyverse)
library(ggprism)

# Load data
data <- read.csv(#path here, check.names = FALSE)

data_filtered <- data %>%
  filter(age %in% c("10w", "20w"))

# Pivot long format
data_long <- data_filtered %>%
  pivot_longer(cols = `0-5pA`:`50-100pA`, 
               names_to = "bin", 
               values_to = "count")

# Set correct bin order
bin_order <- c("0-5pA", "5-10pA", "10-15pA", "15-20pA", "20-25pA", 
               "25-30pA", "30-35pA", "35-40pA", "40-45pA", "45-50pA", 
               "50-100pA", "100+pA")

data_long <- data_long %>%
  mutate(bin = factor(bin, levels = bin_order))

data_long <- data_long %>%
  mutate(group = paste0(genotype, "_", age))

# Colors
custom_colors <- c(
  "wt_10w" = "black",
  "wt_20w" = "black",
  "ko_10w" = "#34d6fa",
  "ko_20w" = "#0213f7"
)

plot_bin_5 <- function(df, age_filter) {
  filtered_df <- df %>%
    filter(age == age_filter) %>%
    mutate(group = factor(paste0(genotype, "_", age), 
                          levels = c("wt_10w", "ko_10w", "wt_20w", "ko_20w")))
  
  ggplot(filtered_df, aes(x = bin, y = count, color = group)) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 8,
      position = position_dodge(width = 0.7)
    ) +
    stat_summary(
      fun.data = function(x) {
        m <- mean(x, na.rm = TRUE)
        se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
        data.frame(y = m, ymin = m - se, ymax = m + se)
      },
      geom = "errorbar",
      width = 0.18,
      linewidth = 0.8,
      position = position_dodge(width = 0.7)
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.7),
      size = 1.2,
      shape = 16
    ) +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0))) +
    labs(
      title = paste("Event counts per bin at", age_filter),
      x = "Current Bin",
      y = "Mean Event Count"
    ) +
    theme_prism(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      legend.position = "none"
    )
}

# Plots
plot_bin_5_10w <- plot_bin_5(data_long, "10w")
plot_bin_5_20w <- plot_bin_5(data_long, "20w")

# Show plots
print(plot_bin_5_usual_10w)
print(plot_bin_5_usual_20w)

# Define output directory
output_dir <- #path

# Make sure directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save function
save_plot <- function(plot_full, filename_base) {
  
  # EPS version
  plot_clean <- plot_full +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      plot.title = element_blank(),
      text = element_text(family = "Arial")  # fix EPS font issue
    )
  
    # Save EPS
  ggsave(filename = file.path(output_dir, paste0(filename_base, "_crossbar.eps")),
         plot = plot_clean,
         device = "eps",
         width = 8,
         height = 6)
  
  # PNG version
  ggsave(filename = file.path(output_dir, paste0(filename_base, "_crossbar.png")),
         plot = plot_full,
         device = "png",
         width = 8,
         height = 6,
         dpi = 300)
}

# Save plots
save_plot(plot_bin_5_10w, "plot_bin_5_10w")
save_plot(plot_bin_5_20w, "plot_bin_5_20w")
