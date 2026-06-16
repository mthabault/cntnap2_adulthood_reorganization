library(readr)
library(dplyr)
library(mgcv)


# Load data
file_sholl <- #path

sholl_data <- read_csv(file_sholl, show_col_types = FALSE) %>%
  mutate(
    mouse_id      = factor(mouse_id),
    cell_id       = factor(cell_id),
    genotype      = factor(genotype),
    age           = factor(age),
    radius        = as.numeric(radius),
    intersections = as.numeric(intersections),
    cell_uid      = interaction(mouse_id, cell_id, drop = TRUE),
    group         = interaction(genotype, age, drop = TRUE)
  ) %>%
  filter(
    genotype %in% c("wt","ko"),
    !is.na(radius),
    !is.na(intersections)
  )

cat("\n========================\n")
cat("Dataset loaded\n")
cat("========================\n")
cat("Rows:", nrow(sholl_data), "\n")
cat("Cells:", nlevels(sholl_data$cell_uid), "\n")


cat("\n========================================\n")
cat("Testing 4-group Sholl smooth differences\n")
cat("========================================\n")

# shared smooth
m0 <- gam(
  intersections ~ genotype * age +
    s(radius, k = 10) +
    s(cell_uid, bs="re"),
  data = sholl_data,
  family = nb(),
  method = "ML"
)

# group-specific smooths
m1 <- gam(
  intersections ~ genotype * age +
    s(radius, by = group, k = 10) +
    s(cell_uid, bs="re"),
  data = sholl_data,
  family = nb(),
  method = "ML"
)

cat("\n--- Model comparison (Chi-square test) ---\n")
cmp <- anova(m0, m1, test="Chisq")
print(cmp)

cat("\n--- Summary of model ---\n")
s <- summary(m1)
print(s)

cat("\n--- Smooth table ---\n")
print(s$s.table)

cat("\nFinished.\n")