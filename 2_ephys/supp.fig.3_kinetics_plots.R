# Load libraries
library(ggplot2)
library(ggprism)
library(dplyr)
library(Cairo)
library(readr)

# Declare Arial font
windowsFonts(Arial = windowsFont("TT Arial"))

# Load data
ipsc_data <- read_csv(
  "#path"
)

# Create group column
ipsc_data <- ipsc_data %>%
  mutate(group = paste(genotype, age, sep = "_"))

# Set factor order
ipsc_data$group <- factor(
  ipsc_data$group,
  levels = c("wt_10w", "ko_10w", "wt_20w", "ko_20w")
)

# Define colors
group_colors <- c(
  "wt_10w" = "grey",
  "ko_10w" = "#34d6fa",
  "wt_20w" = "black",
  "ko_20w" = "#0213f7"
)

# Output path
output_path <- "#path"

# Variables and y-axis limits
plot_info <- list(
  rise = c(0, 10),
  decay = c(0, 40),
  halfwidth = c(0, 30),
  amplitude = c(0, 80)
)

# Plotting function
make_plot <- function(variable, y_limits) {
  
  plot_data <- ipsc_data %>%
    mutate(value = .data[[variable]])
  
  summary_data <- plot_data %>%
    group_by(group) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sem = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups = "drop"
    )
  
  plot_full <- ggplot(plot_data, aes(x = group, y = value)) +
    
    geom_jitter(
      aes(fill = group),
      width = 0.12,
      shape = 21,
      size = 2,
      color = "black",
      alpha = 1
    ) +
    
    geom_crossbar(
      data = summary_data,
      aes(
        x = group,
        y = mean,
        ymin = mean,
        ymax = mean,
        color = group
      ),
      width = 0.35,
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    
    geom_errorbar(
      data = summary_data,
      aes(
        x = group,
        ymin = mean - sem,
        ymax = mean + sem,
        color = group
      ),
      width = 0.15,
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    
    theme_prism(base_size = 14) +
    theme(
      text = element_text(family = "Arial"),
      legend.position = "none"
    ) +
    labs(
      title = paste(variable, "comparison across groups"),
      x = "Group",
      y = variable
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_y_continuous(
      limits = y_limits,
      expand = c(0, 0)
    )
  
  plot_clean <- plot_full +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(
      axis.title = element_blank(),
      legend.position = "none",
      plot.title = element_blank()
    )
  
  print(plot_clean)
  print(plot_full)
  
  ggsave(
    filename = file.path(output_path, paste0(variable, "_clean.eps")),
    plot = plot_clean,
    device = cairo_ps,
    width = 5,
    height = 5,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(output_path, paste0(variable, "_full.png")),
    plot = plot_full,
    width = 5,
    height = 5,
    dpi = 300
  )
  
  cat(variable, "plots saved.\n")
}

# Run all plots
for (var in names(plot_info)) {
  make_plot(var, plot_info[[var]])
}