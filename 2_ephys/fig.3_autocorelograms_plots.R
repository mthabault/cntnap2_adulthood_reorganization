library(tidyverse)
library(ggprism)

output_dir <- #paste path here
data <- read_csv(#paste path here)

# Arrange data & compute mean + SEM
data <- data %>%
  mutate(condition = paste(genotype, age))
lag_levels <- sort(unique(data$lag_ms))

data <- data %>%
  mutate(lag_ms = factor(lag_ms, levels = lag_levels))

summary_data <- data %>%
  group_by(condition, lag_ms) %>%
  summarise(Mean_Autocorr = mean(autocorr_value, na.rm = TRUE),
            SE = sd(autocorr_value, na.rm = TRUE) / sqrt(n()), .groups = "drop")
summary_data <- summary_data %>%
  filter(as.numeric(as.character(lag_ms)) <= 500)

summary_data <- summary_data %>%
  mutate(condition = factor(condition, levels = c("wt 10w", "ko 10w", "wt 20w", "ko 20w")))


# Colors
custom_colors <- c(
  "wt 10w" = "black",
  "wt 20w" = "black",
  "ko 10w" = "#34d6fa",
  "ko 20w" = "#0213f7"
)

# Plotting function
plot_autocorr <- function(data, conditions, title = "", clean = FALSE, custom_colors) {
  
  p <- ggplot(filter(data, condition %in% conditions),
              aes(x = as.numeric(as.character(lag_ms)), y = Mean_Autocorr, color = condition, fill = condition)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = Mean_Autocorr - SE, ymax = Mean_Autocorr + SE), alpha = 0.2, colour = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_x_continuous(limits = c(0, 500), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(-0.01, 0.10), expand = c(0, 0)) +
    theme_prism(base_size = 14) +
    theme(legend.position = "right")
  
  if (!clean) {
    p <- p + labs(title = title, x = "Lag (ms)", y = "Normalized Count (AU)")
  } else {
    p <- p + labs(title = NULL, x = NULL, y = NULL)
  }
  
  return(p)
}


plot_ko_wt_10w_clean <- plot_autocorr(summary_data, c("wt 10w", "ko 10w"), "Autocorrelogram - KO vs WT 10w", clean = TRUE, custom_colors)
plot_ko_wt_10w_full <- plot_autocorr(summary_data, c("wt 10w", "ko 10w"), "Autocorrelogram - KO vs WT 10w", clean = FALSE, custom_colors)

plot_ko_wt_20w_clean <- plot_autocorr(summary_data, c("wt 20w", "ko 20w"), "Autocorrelogram - KO vs WT 20w", clean = TRUE, custom_colors)
plot_ko_wt_20w_full <- plot_autocorr(summary_data, c("wt 20w", "ko 20w"), "Autocorrelogram - KO vs WT 20w", clean = FALSE, custom_colors)

plot_ko_age_clean <- plot_autocorr(summary_data, c("ko 10w", "ko 20w"), "Autocorrelogram - KO 10w vs 20w", clean = TRUE, custom_colors)
plot_ko_age_full <- plot_autocorr(summary_data, c("ko 10w", "ko 20w"), "Autocorrelogram - KO 10w vs 20w", clean = FALSE, custom_colors)


print(plot_ko_wt_10w_full)
print(plot_ko_wt_20w_full)
print(plot_ko_age_full)


ggsave(file.path(output_dir, "autocorr_ko_wt_10w_clean.eps"), plot_ko_wt_10w_clean, width = 8, height = 5, device = cairo_ps, fallback_resolution = 600)
ggsave(file.path(output_dir, "autocorr_ko_wt_10w_full.png"), plot_ko_wt_10w_full, width = 8, height = 5, dpi = 600)

ggsave(file.path(output_dir, "autocorr_ko_wt_20w_clean.eps"), plot_ko_wt_20w_clean, width = 8, height = 5, device = cairo_ps, fallback_resolution = 600)
ggsave(file.path(output_dir, "autocorr_ko_wt_20w_full.png"), plot_ko_wt_20w_full, width = 8, height = 5, dpi = 600)

ggsave(file.path(output_dir, "autocorr_ko_age_clean.eps"), plot_ko_age_clean, width = 8, height = 5, device = cairo_ps, fallback_resolution = 600)
ggsave(file.path(output_dir, "autocorr_ko_age_full.png"), plot_ko_age_full, width = 8, height = 5, dpi = 600)


custom_colors_wt <- c(
  "wt 10w" = "gray60",
  "wt 20w" = "black"
)

plot_wt_age_clean <- plot_autocorr(summary_data, c("wt 10w", "wt 20w"), "Autocorrelogram - WT 10w vs 20w", clean = TRUE, custom_colors_wt)
plot_wt_age_full <- plot_autocorr(summary_data, c("wt 10w", "wt 20w"), "Autocorrelogram - WT 10w vs 20w", clean = FALSE, custom_colors_wt)

print(plot_wt_age_full)

ggsave(file.path(output_dir, "autocorr_wt_age_clean.eps"), plot_wt_age_clean, width = 8, height = 5, device = cairo_ps, fallback_resolution = 600)
ggsave(file.path(output_dir, "autocorr_wt_age_full.png"), plot_wt_age_full, width = 8, height = 5, dpi = 600)

cat("All plots saved.\n")
