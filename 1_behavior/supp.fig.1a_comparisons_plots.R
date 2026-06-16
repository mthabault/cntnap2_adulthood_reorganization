library(dplyr)
library(ggplot2)
library(readr)
library(ggprism)
library(Cairo)

behaviour_data <- read.csv(#paste path here,
  stringsAsFactors = FALSE
)

#Create composite behaviors: repetitive & stereotyped / exploratory
composites <- behaviour_data %>%
  mutate(
    composite = case_when(
      behavior %in% c("scratching", "head_body_twitch", "grooming_long", "grooming_short") ~
        "repetitive and stereotyped behaviour",
      behavior %in% c("digging", "supported_rearing", "unsupported_rearing") ~
        "exploratory behaviour",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(composite)) %>%
  group_by(video, mouse, genotype, age, composite) %>%
  summarise(
    n_behaviours = sum(n_behaviours, na.rm = TRUE),
    total_duration = sum(total_duration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(behavior = composite)

behaviour_data2 <- bind_rows(behaviour_data, composites)

behaviours <- c(
  "repetitive and stereotyped behaviour",
  "exploratory behaviour"
)

behaviour_data2 <- behaviour_data2 %>%
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

#define colors
point_colors <- c(
  "wt_10w" = "white",
  "ko_10w" = "#34d6fa",
  "wt_20w" = "grey",
  "ko_20w" = "#0213f7"
)

#plots
make_spaghetti_plot <- function(data, behaviour_name) {
  
  df <- data %>%
    dplyr::filter(behavior == behaviour_name)
  
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
    scale_x_discrete(expand = expansion(add = 1.5)) +
    scale_y_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, by = 20),
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

for (b in behaviours) {
  cat("\n\n", b, "\n")
  p <- make_spaghetti_plot(behaviour_data2, b)
  print(p)
}

# Generate and save the 7 plots as EPS
output_dir <- #paste path here

for (b in behaviours) {
  cat("\n\n===== ", b, "\n")
  
  p <- make_spaghetti_plot(behaviour_data2, b)
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