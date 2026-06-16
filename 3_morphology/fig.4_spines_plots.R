# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(ggprism)

# Load the datasets
spines <- read_csv(#path)
dendrites <- read_csv(#path)

dendrites_len <- dendrites %>%
  select(mouse_id, cell_id, ID, dendrite_length) %>%
  rename(
    parent_id = ID,
    segment_length = dendrite_length
  )

spines_seg <- spines %>%
  filter(depth %in% c(1, 2, 3, 4)) %>%
  left_join(
    dendrites_len,
    by = c("mouse_id", "cell_id", "parent_id"),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(segment_length))

segment_counts <- spines_seg %>%
  group_by(mouse_id, cell_id, genotype, age, parent_id, depth, segment_length) %>%
  summarise(
    n_total = n(),
    n_thin  = sum(spine_morphology == "Thin", na.rm = TRUE),
    n_stub  = sum(spine_morphology == "Stubby", na.rm = TRUE),
    n_mush  = sum(spine_morphology == "Mushroom", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    density_total = (n_total / segment_length) * 10,
    density_thin  = (n_thin / segment_length) * 10,
    density_stub  = (n_stub / segment_length) * 10,
    density_mush  = (n_mush / segment_length) * 10,
    order         = depth + 1,
    order_label   = factor(order),
    group         = paste(genotype, age, sep = "_")
  )

# Colors
group_colors <- c(
  "wt_10w" = "grey",
  "ko_10w" = "#34d6fa",
  "wt_20w" = "black",
  "ko_20w" = "#0213f7"
)

# Plotting function
plot_one_order_one_density <- function(df, order_value, yvar, ylab) {
  
  df_plot <- df %>%
    filter(order == order_value) %>%
    mutate(group = factor(group, levels = names(group_colors)))
  
  summary_data <- df_plot %>%
    group_by(group) %>%
    summarise(
      mean = mean(.data[[yvar]], na.rm = TRUE),
      sem = sd(.data[[yvar]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[yvar]]))),
      .groups = "drop"
    )
  
  if (yvar == "density_total") {
    y_max <- 15
    y_breaks <- seq(0, 15, 5)
  } else {
    y_max <- 10
    y_breaks <- seq(0, 10, 2)
  }
  
  p <- ggplot(df_plot, aes(x = group, y = .data[[yvar]])) +
    
    geom_jitter(
      aes(fill = group),
      width = 0.12,
      height = 0.08,
      shape = 21,
      size = 2,
      color = "black",
      alpha = 1
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
    
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(ylim = c(0, y_max)) +
    
    labs(
      x = NULL,
      y = ylab,
      title = paste("Order", order_value)
    ) +
    theme_prism() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  p
}

# Density labels
densities <- list(
  density_total = "Total spine density (count / 10 µm)",
  density_mush  = "Mushroom spine density (count / 10 µm)",
  density_stub  = "Stubby spine density (count / 10 µm)",
  density_thin  = "Thin spine density (count / 10 µm)"
)

# Generate plots
plots <- list()

for (o in 2:5) {
  for (nm in names(densities)) {
    key <- paste0(nm, "_order_", o)
    plots[[key]] <- plot_one_order_one_density(
      segment_counts,
      order_value = o,
      yvar = nm,
      ylab = densities[[nm]]
    )
  }
}

# Print all plots
for (p in plots) print(p)

output_dir <- #path

for (name in names(plots)) {
  
  cairo_ps(
    filename = file.path(output_dir, paste0(name, ".eps")),
    width = 4,
    height = 4,
    onefile = FALSE
  )
  
  print(
    plots[[name]] +
      theme(text = element_text(family = "Helvetica"))
  )
  
  dev.off()
}