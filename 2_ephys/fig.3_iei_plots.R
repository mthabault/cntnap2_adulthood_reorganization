# Load required libraries
library(ggplot2)
library(ggprism)
library(tidyr)
library(dplyr)
library(Cairo)

# Declare Arial font properly for plotting
windowsFonts(Arial = windowsFont("TT Arial"))

# Load the data
iei_data <- read.csv(#path here)

# Reshape from wide to long format
iei_long <- iei_data %>%
  pivot_longer(cols = starts_with("X"), names_to = "iei_bin", values_to = "count") %>%
  mutate(iei_bin = as.numeric(gsub("X", "", iei_bin)))  # Convert "X100" to numeric 100

# Define plotting function
plot_iei <- function(data, target_age, colors, title_suffix, clean = FALSE) {
  
  # Filter for target age and genotypes
  filtered <- data %>%
    filter(age == target_age, genotype %in% c("wt", "ko"))
  
  # Summarise mean and SEM per bin
  summary <- filtered %>%
    group_by(genotype, age, iei_bin) %>%
    summarise(
      mean_count = mean(count, na.rm = TRUE),
      sem_count = sd(count, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Build the plot
  p <- ggplot(summary, aes(x = iei_bin, y = mean_count, color = interaction(genotype, age), group = interaction(genotype, age))) +
    geom_point(size = 2) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = mean_count - sem_count, ymax = mean_count + sem_count), width = 100, linewidth = 0.8) +
    theme_prism() +
    theme(text = element_text(family = "Arial")) +  # Compatible font
    scale_x_continuous(expand = c(0, 0), limits = c(0, 5000)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_color_manual(values = colors)
  
  # Add or remove labels and legend based on "clean" argument
  if (!clean) {
    p <- p + labs(
      title = paste("IEI Mean ± SEM: WT vs KO (", title_suffix, ")", sep = ""),
      x = "Inter-Event Interval Bin (ms)",
      y = "Mean Count",
      color = "Genotype x Age"
    )
  } else {
    p <- p + labs(title = NULL, x = NULL, y = NULL, color = NULL) +
      theme(legend.position = "none")
  }
  
  return(p)
}

# Generate plots for 10w and 20w
p10_clean <- plot_iei(iei_long, "10w", c("wt.10w" = "black", "ko.10w" = "#34d6fa"), "10w", clean = TRUE)
p10_full <- plot_iei(iei_long, "10w", c("wt.10w" = "black", "ko.10w" = "#34d6fa"), "10w", clean = FALSE)

p20_clean <- plot_iei(iei_long, "20w", c("wt.20w" = "black", "ko.20w" = "#0213f7"), "20w", clean = TRUE)
p20_full <- plot_iei(iei_long, "20w", c("wt.20w" = "black", "ko.20w" = "#0213f7"), "20w", clean = FALSE)

print(p10_full)
print(p20_full)

# Define output directory
output_dir <- #path

# Save clean version as EPS
ggsave(filename = file.path(output_dir, "IEI_10w_clean.eps"),
       plot = p10_clean,
       width = 8, height = 5, device = cairo_ps, fallback_resolution = 600)

ggsave(filename = file.path(output_dir, "IEI_20w_clean.eps"),
       plot = p20_clean,
       width = 8, height = 5, device = cairo_ps, fallback_resolution = 600)

# Save full version as PNG
ggsave(filename = file.path(output_dir, "IEI_10w_full.png"),
       plot = p10_full,
       width = 8, height = 5, dpi = 600)

ggsave(filename = file.path(output_dir, "IEI_20w_full.png"),
       plot = p20_full,
       width = 8, height = 5, dpi = 600)

cat("All EPS and PNG files saved to:\n", output_dir, "\n")
