# Load libraries
library(ggplot2)
library(ggprism)
library(dplyr)
library(Cairo)

# Declare Arial for Windows
windowsFonts(Arial = windowsFont("TT Arial"))

behaviour_data <- read.csv(#paste path here)

# define list of behaviors
behaviours <- c(
  "digging",
  "supported_rearing",
  "unsupported_rearing",
  "scratching",
  "head_body_twitch",
  "grooming_short",
  "grooming_long"
)

# define colors
point_colors <- c(
  "wt_10w" = "white",
  "ko_10w" = "#34d6fa",
  "wt_20w" = "grey",
  "ko_20w" = "#0213f7"
)

# Function for specific y axis (plotting preference, can be removed)
compute_yaxis <- function(y) {
  max_val <- max(y, na.rm = TRUE)
  ymax <- ceiling(max_val / 10) * 10 + 10
  if (ymax == 0) ymax <- 10
  
  by_step <- if (ymax <= 10) {
    5
  } else if (ymax <= 50) {
    10
  } else if (ymax <= 100) {
    20
  } else if (ymax <= 200) {
    50
  } else {
    50
  }
  
  list(ymax = ymax, by_step = by_step)
}

# Clean and prepare data
behaviour_data <- behaviour_data %>%
  filter(
    behavior %in% behaviours,
    age %in% c("10w", "20w")
  ) %>%
  mutate(
    genotype = tolower(genotype),
    genotype = factor(genotype, levels = c("wt", "ko")),
    age      = factor(age, levels = c("10w", "20w")),
    group    = factor(
      paste(genotype, age, sep = "_"),
      levels = c("wt_10w", "ko_10w", "wt_20w", "ko_20w")
    )
  )

# Plots
make_spaghetti_plot <- function(data, behaviour_name) {
  
  df <- data %>%
    dplyr::filter(behavior == behaviour_name)
  
  # y-axis limits (manual per behaviour)
  if (behaviour_name %in% c("scratching", "head_body_twitch",
                            "grooming_long", "grooming_short")) {
    yl <- list(ymax = 30, by_step = 10)
  } else if (behaviour_name %in% c("digging", "supported_rearing",
                                   "unsupported_rearing")) {
    yl <- list(ymax = 60, by_step = 10)
  } else {
    yl <- compute_yaxis(df$n_behaviours)  # fallback (should never be used)
  }
  
  
  # build label strings: one row per unique (genotype, age, n_behaviours)
  df_labels <- df %>%
    dplyr::group_by(genotype, age, n_behaviours) %>%
    dplyr::summarise(
      mouse_label = paste(sort(unique(mouse)), collapse = "; "),
      .groups = "drop"
    )
  
  ggplot(df, aes(x = age, y = n_behaviours, group = mouse)) +
    geom_line(color = "black", alpha = 1, linewidth = 1.5) +
    geom_point(
      aes(fill = group),
      size  = 2.5,
      shape = 21,
      color = "black",
      stroke = 0.6
    ) +
    
    geom_text(
      data = df_labels %>% dplyr::filter(age == "10w"),
      aes(x = age, y = n_behaviours, label = mouse_label),
      inherit.aes = FALSE,
      hjust  = 1,
      nudge_x = -0.08,
      vjust  = 0.5,
      size   = 3,
      family = "Arial"
    ) +
    
    geom_text(
      data = df_labels %>% dplyr::filter(age == "20w"),
      aes(x = age, y = n_behaviours, label = mouse_label),
      inherit.aes = FALSE,
      hjust  = 0,
      nudge_x = 0.08,
      vjust  = 0.5,
      size   = 3,
      family = "Arial"
    ) +
    
    scale_fill_manual(values = point_colors) +
    
    scale_x_discrete(
      expand = expansion(add = 1.5)
    ) +
    
    scale_y_continuous(
      limits = c(0, yl$ymax),
      breaks = seq(0, yl$ymax, by = yl$by_step),
      expand = expansion(mult = c(0.05, 0.05)) 
    ) +
    
    facet_wrap(
      ~ genotype,
      nrow = 1,
      labeller = as_labeller(c("wt" = "WT", "ko" = "KO"))
    ) +
    
    labs(
      title = behaviour_name,
      x = "",
      y = "Number of events"
    ) +
    
    theme_prism(base_size = 16) +
    theme(
      text        = element_text(family = "Arial"),
      plot.title  = element_text(hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(10, 60, 10, 60)
    )
}

# Generate and print the 7 plots
for (b in behaviours) {
  cat("\n\n===== ", b, " – n_behaviours (10w vs 20w, WT/KO) =====\n")
  p <- make_spaghetti_plot(behaviour_data, b)
  print(p)
}

# Generate and save the 7 plots as EPS
output_dir <- #paste path here

for (b in behaviours) {
  cat("\n\n===== ", b, " – n_behaviours (10w vs 20w, WT/KO) =====\n")
  
  p <- make_spaghetti_plot(behaviour_data, b)
  print(p)
  
  # Save each plot as EPS
  eps_filename <- file.path(output_dir, paste0(b, "_spaghetti.eps"))
  ggsave(
    filename = eps_filename,
    plot = p,
    device = "eps",
    width = 8,
    height = 4,
    dpi = 300
  )
}


