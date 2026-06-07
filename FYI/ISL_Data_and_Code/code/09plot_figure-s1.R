# Load required R packages for data manipulation, scientific visualization, and statistical annotation
library(tidyverse)    # Core package for data cleaning (dplyr) and visualization (ggplot2)
library(ggsci)        # For scientific journal-style color palettes (not explicitly used here but loaded)
library(paletteer)    # For accessing additional color palettes (not explicitly used here but loaded)
library(cartography)  # For spatial data visualization (not used in this script, likely a leftover)
library(ggpubr)       # For enhancing publication-quality plots (e.g., statistical tests)
library(ggpmisc)      # For adding polynomial regression equations, R², and p-values to plots

# --- Data Loading and Preprocessing ---
# Read community data (contains anuran richness/abundance and biogeographic predictors)
# Apply log10 transformations to skewed biogeographic variables to meet linear model assumptions
dat <- read_csv("data/community_data.csv") %>%
  mutate(
    isolation = log10(isolation / 1000),    # Log10-transform isolation (convert to km first)
    area = log10(area),                     # Log10-transform island area
    water_area = log10(water_area)          # Log10-transform water body area
  )

# Convert 'invasion' (invasion status: "N" = non-invaded, "Y" = invaded) to factor for grouping
dat$invasion <- as.factor(dat$invasion)

# Create data subsets for comparing Control (non-invaded) vs. Invaded groups:
# 1. Control group: Only non-invaded ponds (invasion == "N"), labeled as "Control"
dt_control <- dat %>% filter(invasion == "N") %>% mutate(type = "Control")
# 2. Invaded group: All ponds labeled as "Invaded" (for group-wise fitting)
dt_addinv <- dat %>% mutate(type = "Invaded")
# 3. Combined dataset: Merge Control and Invaded subsets for unified model fitting
dt_all <- rbind(dt_control, dt_addinv)
# Set factor levels for 'type' to ensure "Control" appears before "Invaded" in plots/legends
dt_all$type <- factor(dt_all$type, levels = c("Control", "Invaded"))


# --- Section 1: Native Anuran Richness vs. Biogeographic Predictors ---
# Goal: Explore how native frog species richness relates to isolation, island area, and water area
# (Compare patterns between Control and Invaded groups)

# figs1a: Isolation (log10 km) vs. Native Anuran Richness
figs1a <- ggplot() +
  # Scatter points: All ponds labeled as "Invaded" (dt_addinv), colored by actual invasion status
  geom_point(
    data = dt_addinv,
    aes(x = isolation, y = native_frog_richness, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10  # Shape 21 = filled circle with border; alpha = transparency (50%)
  ) +
  # Axis labels: Y = Native anuran richness (note: typo in original code: "richnes" → corrected in label), X = Isolation
  labs(y = "Native anuran richness", x = "Distance to closest mainland [log10(km)]") +
  # Color/fill scale: Blue (#3498DB) for non-invaded ("N"), Orange (#FF8000) for invaded ("Y")
  # 4 values provided (redundant, only 2 used) for consistency with other plots
  scale_fill_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  # Group-wise linear regression fits (method = "lm") for Control vs. Invaded
  stat_smooth(
    data = dt_all,
    aes(x = isolation, y = native_frog_richness, color = type, fill = type),
    method = "lm", size = 1.5, linetype = 2  # Linetype 2 = dashed line
  ) +
  # Add regression equation, R², and p-value (formatted with "~~~~" for line breaks in plotmath)
  stat_poly_eq(
    data = dt_all,
    aes(
      x = isolation, y = native_frog_richness, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")
    ),
    size = 3  # Text size for annotations
  )


# figs1b: Island Area (log10 km²) vs. Native Anuran Richness
figs1b <- ggplot() +
  geom_point(
    data = dt_addinv,
    aes(x = area, y = native_frog_richness, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Native anuran richness", x = "Island area [log10(km2)]") +
  scale_fill_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  stat_smooth(
    data = dt_all,
    aes(x = area, y = native_frog_richness, color = type, fill = type),
    method = "lm", size = 1.5, linetype = 2
  ) +
  stat_poly_eq(
    data = dt_all,
    aes(
      x = area, y = native_frog_richness, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")
    ),
    size = 3
  )


# figs1c: Water Body Area (log10 km²) vs. Native Anuran Richness
figs1c <- ggplot() +
  geom_point(
    data = dt_addinv,
    aes(x = water_area, y = native_frog_richness, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Native anuran richness", x = "Water body area [log10(km2)]") +
  scale_fill_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  stat_smooth(
    data = dt_all,
    aes(x = water_area, y = native_frog_richness, color = type, fill = type),
    method = "lm", size = 1.5, linetype = 2
  ) +
  stat_poly_eq(
    data = dt_all,
    aes(
      x = water_area, y = native_frog_richness, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")
    ),
    size = 3
  )


# --- Section 2: Native Anuran Abundance vs. Biogeographic Predictors ---
# Goal: Explore how native frog abundance relates to isolation, island area, and water area
# (Same structure as Section 1, but response variable = native_frog_abundance)

# figs1d: Isolation (log10 km) vs. Native Anuran Abundance
figs1d <- ggplot() +
  geom_point(
    data = dt_addinv,
    aes(x = isolation, y = native_frog_abundance, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Native anuran abundance", x = "Distance to closest mainland [log10(km)]") +
  scale_fill_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  stat_smooth(
    data = dt_all,
    aes(x = isolation, y = native_frog_abundance, color = type, fill = type),
    method = "lm", size = 1.5, linetype = 2
  ) +
  stat_poly_eq(
    data = dt_all,
    aes(
      x = isolation, y = native_frog_abundance, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")
    ),
    size = 3
  )


# figs1e: Island Area (log10 km²) vs. Native Anuran Abundance
figs1e <- ggplot() +
  geom_point(
    data = dt_addinv,
    aes(x = area, y = native_frog_abundance, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Native anuran abundance", x = "Island area [log10(km2)]") +
  scale_fill_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  stat_smooth(
    data = dt_all,
    aes(x = area, y = native_frog_abundance, color = type, fill = type),
    method = "lm", size = 1.5, linetype = 2
  ) +
  stat_poly_eq(
    data = dt_all,
    aes(
      x = area, y = native_frog_abundance, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")
    ),
    size = 3
  )


# figs1f: Water Body Area (log10 km²) vs. Native Anuran Abundance
figs1f <- ggplot() +
  geom_point(
    data = dt_addinv,
    aes(x = water_area, y = native_frog_abundance, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Native anuran abundance", x = "Water body area [log10(km2)]") +
  scale_fill_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000", "#3498DB", "#FF8000")) +
  stat_smooth(
    data = dt_all,
    aes(x = water_area, y = native_frog_abundance, color = type, fill = type),
    method = "lm", size = 1.5, linetype = 2
  ) +
  stat_poly_eq(
    data = dt_all,
    aes(
      x = water_area, y = native_frog_abundance, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")
    ),
    size = 3
  )


# --- Section 3: Bullfrog Abundance vs. Biogeographic Predictors (Invaded ponds Only) ---
# Goal: Explore how invasive bullfrog abundance relates to biogeographic predictors
# (Only uses invaded ponds: dt_invaded = dat %>% filter(invasion == "Y"))

# Subset data to only invaded ponds
dt_invaded <- dat %>% filter(invasion == "Y")

# figs1g: Isolation (log10 km) vs. Bullfrog Abundance (Invaded ponds Only)
figs1g <- ggplot(data = dt_invaded, aes(x = isolation, y = bullfrog_abundance)) +
  # Scatter points: Uniform orange color (#FF8000) for all invaded ponds
  geom_point(size = 3, shape = 21, alpha = 5/10, fill = "#FF8000", colour = "#FF8000") +
  # Linear regression fit (only for invaded ponds)
  stat_smooth(method = "lm", size = 1.5, linetype = 2, fill = "#FF8000", colour = "#FF8000") +
  # Add regression equation, R², and p-value
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")),
    size = 3
  ) +
  labs(y = "Bullfrog abundance", x = "Distance to closest mainland [log10(km)]")


# figs1h: Island Area (log10 km²) vs. Bullfrog Abundance (Invaded ponds Only)
figs1h <- ggplot(data = dt_invaded, aes(x = area, y = bullfrog_abundance)) +
  geom_point(size = 3, shape = 21, alpha = 5/10, fill = "#FF8000", colour = "#FF8000") +
  stat_smooth(method = "lm", size = 1.5, linetype = 2, fill = "#FF8000", colour = "#FF8000") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")),
    size = 3
  ) +
  labs(y = "Bullfrog abundance", x = "Island area [log10(km2)]")


# figs1i: Water Body Area (log10 km²) vs. Bullfrog Abundance (Invaded ponds Only)
figs1i <- ggplot(data = dt_invaded, aes(x = water_area, y = bullfrog_abundance)) +
  geom_point(size = 3, shape = 21, alpha = 5/10, fill = "#FF8000", colour = "#FF8000") +
  stat_smooth(method = "lm", size = 1.5, linetype = 2, fill = "#FF8000", colour = "#FF8000") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")),
    size = 3
  ) +
  labs(y = "Bullfrog abundance", x = "Water body area [log10(km2)]")


# --- Assemble and Save Combined Figure (figs1) ---
# Combine all 9 subfigures (figs1a-i) into a 3-column, 3-row grid
# labels = c("AUTO"): Automatically generates subfigure labels (A, B, C, ..., I)
figs1 <- cowplot::plot_grid(
  figs1a, figs1b, figs1c, figs1d, figs1e, figs1f, figs1g, figs1h, figs1i,
  ncol = 3,  # 3 columns (rows = 3, since 9 plots total)
  labels = c("AUTO")  # Auto-generate alphabetical labels (A-I)
  # Original comment: Extra labels ("D", "E", "F") are commented out (not needed)
)

# Save the combined figure as a high-resolution PNG (suitable for publications)
# Width = 9 inches, Height = 9 inches (balances readability for 3x3 grid)
ggsave(
  filename = paste0('figure/figure_s1_community_biogeography.png'),
  plot = figs1,
  width = 9,
  height = 9
)
