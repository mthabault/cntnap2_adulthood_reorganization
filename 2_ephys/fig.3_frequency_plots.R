# Load required libraries
library(ggplot2)
library(ggprism)
library(dplyr)
library(Cairo)

# Arial font for eps export
windowsFonts(Arial = windowsFont("TT Arial"))

# Load the data
ipsc_data <- read.csv(#paste path here)

# Create 'group' column
ipsc_data <- ipsc_data %>%
  mutate(group = paste(genotype, age, sep = "_"))

# Define plotting function
plot_frequency <- function(data, target_age, colors, title_suffix) {

  filtered <- data %>%
    filter(age == target_age, genotype %in% c("wt", "ko")) %>%
    mutate(group = factor(paste(genotype, age, sep = "_"),
                          levels = c(paste0("wt_", target_age), paste0("ko_", target_age))))
  
  ggplot(filtered, aes(x = group, y = frequency, fill = group)) +
    stat_summary(
      aes(color = group),
      fun = mean,
      geom = "crossbar",
      width = 0.5,
      linewidth = 1,
      middle.linewidth = 1
    ) +
    stat_summary(
      aes(color = group),
      fun.data = function(x) {
        m <- mean(x, na.rm = TRUE)
        se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
        data.frame(y = m, ymin = m - se, ymax = m + se)
      },
      geom = "errorbar",
      width = 0.15,
      linewidth = 1
    ) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0),
      size = 2,
      alpha = 1,
      shape = 21,
      color = "black"
    ) +
    theme_prism(base_size = 14) +
    theme(
      text = element_text(family = "Arial"),
      legend.position = "none"
    ) +
    labs(
      title = paste("freq"),
      x = "Group",
      y = "Frequency"
    ) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(limits = c(0, 4), expand = c(0, 0))
}

# Define colors
colors_10w <- c("wt_10w" = "black", "ko_10w" = "#34d6fa")
colors_20w <- c("wt_20w" = "black", "ko_20w" = "#0213f7")

# Generate plots
freq10 <- plot_frequency(ipsc_data, "10w", colors_10w, "10w")
freq20 <- plot_frequency(ipsc_data, "20w", colors_20w, "20w")

print(freq10)
print(freq20)

# Define output directory
output_dir <- #paste path here
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


make_clean_version <- function(plot) {
  plot +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(
      axis.title = element_blank(),
      legend.position = "none",
      plot.title = element_blank()
    )
}


# clean version = eps to be arranged in figure
freq10_clean <- make_clean_version(freq10)
freq20_clean <- make_clean_version(freq20)

# Save clean as EPS
ggsave(filename = file.path(output_dir, "crossbar_iPSC_Frequency_10w_clean.eps"),
       plot = freq10_clean, width = 3, height = 5, device = cairo_ps, dpi = 600)
ggsave(filename = file.path(output_dir, "crossbar_iPSC_Frequency_20w_clean.eps"),
       plot = freq20_clean, width = 3, height = 5, device = cairo_ps, dpi = 600)

# Save full version as PNG
ggsave(filename = file.path(output_dir, "crossbar_iPSC_Frequency_10w_full.png"),
       plot = freq10, width = 3, height = 5, dpi = 300)
ggsave(filename = file.path(output_dir, "crossbar_iPSC_Frequency_20w_full.png"),
       plot = freq20, width = 3, height = 5, dpi = 300)

Cat("All plots saved to:\n", output_dir, "\n")
