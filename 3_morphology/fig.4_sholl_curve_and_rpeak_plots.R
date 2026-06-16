library(dplyr)
library(readr)
library(ggplot2)
library(ggprism)

file_sholl <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/3_morphology/3_csv_datasets/sholl.csv"

file_metrics <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/0_ready_to_submit/3_morphology/3_csv_datasets/sholl_metrics_per_cell.csv"

output_dir <- "C:/Users/Mathieu_Thabault/Desktop/papier_cnt/dev"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sholl_data <- read_csv(file_sholl) %>%
  mutate(
    genotype  = as.factor(genotype),
    age       = as.factor(age),
    sample_id = paste(mouse_id, cell_id, sep = "_")
  )

metrics <- read_csv(file_metrics) %>%
  mutate(
    genotype = factor(genotype),
    age      = factor(age),
    group    = factor(paste(genotype, age, sep = "_"))
  )

p_to_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")
           )
    )
  )
}

# Sholl curves
plot_sholl <- function(data_wt, data_ko, age_label, color_ko, output_dir) {
  
  data_wt <- data_wt %>% mutate(sample_id = paste(mouse_id, cell_id, sep = "_"))
  data_ko <- data_ko %>% mutate(sample_id = paste(mouse_id, cell_id, sep = "_"))
  
  summary_wt <- data_wt %>%
    group_by(radius) %>%
    summarize(
      mean_intersections = mean(intersections),
      sem_intersections  = sd(intersections) / sqrt(n()),
      .groups = "drop"
    )
  
  summary_ko <- data_ko %>%
    group_by(radius) %>%
    summarize(
      mean_intersections = mean(intersections),
      sem_intersections  = sd(intersections) / sqrt(n()),
      .groups = "drop"
    )
  
  summary_sub <- summary_wt %>%
    inner_join(summary_ko,
               by = "radius",
               suffix = c("_wt", "_ko")) %>%
    mutate(diff_mean = mean_intersections_wt - mean_intersections_ko)
  
  p <- ggplot() +
    geom_ribbon(
      data = summary_sub,
      aes(
        x = radius,
        ymin = pmin(diff_mean, 0),
        ymax = pmax(diff_mean, 0)
      ),
      fill = "#fcaa12",
      alpha = 1,
      colour = NA
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey60",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    { if (age_label == "20w")
      geom_vline(
        xintercept = c(100, 120),
        color = "grey60",
        linewidth = 0.7,
        linetype = "dotted"
      ) else NULL } +
    geom_ribbon(
      data = summary_wt,
      aes(
        x = radius,
        ymin = pmax(mean_intersections - sem_intersections, 0),
        ymax = mean_intersections + sem_intersections
      ),
      fill = "grey50",
      alpha = 1
    ) +
    geom_line(
      data = summary_wt,
      aes(x = radius, y = mean_intersections),
      color = "black",
      linewidth = 1.3
    ) +
    geom_ribbon(
      data = summary_ko,
      aes(
        x = radius,
        ymin = pmax(mean_intersections - sem_intersections, 0),
        ymax = mean_intersections + sem_intersections
      ),
      fill = color_ko,
      alpha = 1
    ) +
    geom_line(
      data = summary_ko,
      aes(x = radius, y = mean_intersections),
      color = color_ko,
      linewidth = 1.3
    ) +
    labs(
      title = paste0("Sholl ", age_label),
      x = "Distance from soma (µm)",
      y = "Intersections"
    ) +
    theme_prism(base_size = 12) +
    scale_y_continuous(
      breaks = seq(-15, 30, by = 5),
      limits = c(-15, 30),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      limits = c(0, 200),
      breaks = seq(0, 200, 50),
      expand = c(0, 0)
    ) +
    theme(text = element_text(family = "Helvetica"))
  
  print(p)
  
  ggsave(
    filename = file.path(output_dir, paste0("sholl_", age_label, "_plot.png")),
    plot = p,
    width = 6, height = 4, units = "in", dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0("sholl_", age_label, "_plot.eps")),
    plot = p,
    width = 6, height = 4, units = "in",
    device = "eps",
    family = "Helvetica"
  )
}

data_wt_10w <- subset(sholl_data, genotype == "wt" & age == "10w")
data_ko_10w <- subset(sholl_data, genotype == "ko" & age == "10w")
plot_sholl(data_wt_10w, data_ko_10w, "10w", "#34d6fa", output_dir)

data_wt_20w <- subset(sholl_data, genotype == "wt" & age == "20w")
data_ko_20w <- subset(sholl_data, genotype == "ko" & age == "20w")
plot_sholl(data_wt_20w, data_ko_20w, "20w", "#0213f7", output_dir)

cat("\nAll Sholl plots saved as PNG + EPS.\n")

ages <- c("10w", "20w")

for (a in ages) {
  cat("\nAge:", a, "\n")
  
  sub   <- metrics %>% filter(age == a)
  t_fit <- t.test(sub$Rmax ~ sub$genotype, var.equal = FALSE)
  
  t_val <- unname(t_fit$statistic)   # t value
  df_t  <- unname(t_fit$parameter)   # Welch df
  p_t   <- t_fit$p.value
  star_t <- p_to_stars(p_t)
  
  cat(
    "Welch t-test for Rmax, WT vs KO:",
    "t =", round(t_val, 3),
    "df =", round(df_t, 1),
    "p =", signif(p_t, 3),
    star_t, "\n"
  )
}


cat("\nDone with Welch tests.\n")

# Colors
colors <- c(
  "wt_10w" = "grey",
  "ko_10w" = "#34d6fa",
  "wt_20w" = "black",
  "ko_20w" = "#0213f7"
)

# Rpeak (= Rmax here) plots
rmax_10w <- metrics %>%
  filter(age == "10w") %>%
  mutate(group_label = factor(paste(genotype, age, sep = "_"),
                              levels = c("wt_10w", "ko_10w")))

rmax_10w_summary <- rmax_10w %>%
  group_by(group_label) %>%
  summarise(
    mean = mean(Rmax, na.rm = TRUE),
    sem = sd(Rmax, na.rm = TRUE) / sqrt(sum(!is.na(Rmax))),
    .groups = "drop"
  )

p_rmax_10w <- ggplot(rmax_10w, aes(x = group_label, y = Rmax)) +
  geom_jitter(
    aes(fill = group_label),
    width = 0.12,
    shape = 21,
    size = 1,
    color = "black",
    alpha = 1
  ) +
  geom_errorbar(
    data = rmax_10w_summary,
    aes(x = group_label, ymin = mean - sem, ymax = mean + sem, color = group_label),
    width = 0.15,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_crossbar(
    data = rmax_10w_summary,
    aes(x = group_label, y = mean, ymin = mean, ymax = mean, color = group_label),
    width = 0.35,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(x = NULL, y = "Rmax (µm)") +
  scale_y_continuous(
    limits = c(0, 200),
    breaks = seq(0, 200, 50),
    expand = c(0, 0)
  ) +
  theme_prism(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

print(p_rmax_10w)

ggsave(
  filename = file.path(output_dir, "rmax_10w_jitter.png"),
  plot = p_rmax_10w,
  width = 3,
  height = 4,
  units = "in",
  dpi = 300
)


rmax_20w <- metrics %>%
  filter(age == "20w") %>%
  mutate(group_label = factor(paste(genotype, age, sep = "_"),
                              levels = c("wt_20w", "ko_20w")))

rmax_20w_summary <- rmax_20w %>%
  group_by(group_label) %>%
  summarise(
    mean = mean(Rmax, na.rm = TRUE),
    sem = sd(Rmax, na.rm = TRUE) / sqrt(sum(!is.na(Rmax))),
    .groups = "drop"
  )

p_rmax_20w <- ggplot(rmax_20w, aes(x = group_label, y = Rmax)) +
  geom_jitter(
    aes(fill = group_label),
    width = 0.12,
    shape = 21,
    size = 1,
    color = "black",
    alpha = 1
  ) +
  geom_errorbar(
    data = rmax_20w_summary,
    aes(x = group_label, ymin = mean - sem, ymax = mean + sem, color = group_label),
    width = 0.15,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_crossbar(
    data = rmax_20w_summary,
    aes(x = group_label, y = mean, ymin = mean, ymax = mean, color = group_label),
    width = 0.35,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(x = NULL, y = "Rmax (µm)") +
  scale_y_continuous(
    limits = c(0, 200),
    breaks = seq(0, 200, 50),
    expand = c(0, 0)
  ) +
  theme_prism(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

print(p_rmax_20w)

ggsave(
  filename = file.path(output_dir, "rmax_20w_jitter.png"),
  plot = p_rmax_20w,
  width = 3,
  height = 4,
  units = "in",
  dpi = 300
)

cat("\nSaving Rmax jitter plots as EPS...\n")

cairo_ps(
  filename = file.path(output_dir, "rmax_10w_jitter.eps"),
  width = 3,
  height = 4,
  onefile = FALSE,
  family = "Helvetica"
)
print(p_rmax_10w)
dev.off()

cairo_ps(
  filename = file.path(output_dir, "rmax_20w_jitter.eps"),
  width = 3,
  height = 4,
  onefile = FALSE,
  family = "Helvetica"
)
print(p_rmax_20w)
dev.off()

cat("\nRmax jitter plots saved as PNG + EPS.\n")