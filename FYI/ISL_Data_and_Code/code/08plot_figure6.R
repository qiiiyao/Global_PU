# Load required R packages for data manipulation, visualization, and plot assembly
library(tidyverse)    # For data cleaning (dplyr) and visualization (ggplot2)
library(cowplot)      # For combining multiple ggplot2 plots into a single cohesive figure
library(ggpubr)       # For adding statistical annotations (e.g., regression equations) to plots

# --- Set Global ggplot2 Theme for Consistent Publication-Ready Formatting ---
# Standardize plot style across all subfigures to ensure uniformity in appearance
theme_set(theme_test() +
            theme(
              panel.grid = element_blank(),          # Remove background grid lines to reduce clutter
              legend.position = "none",              # Hide legend by default (no redundant legends needed)
              legend.title = element_blank(),        # Remove legend title (irrelevant if legend is hidden)
              legend.text = element_text(color = "black", size = 10),  # Legend text style (fallback)
              axis.text.x = element_text(size = 10, color = "black"),   # X-axis tick label style (readability)
              axis.text.y = element_text(size = 10, color = "black"),   # Y-axis tick label style (readability)
              axis.title.y = element_text(size = 10, color = "black", angle = 90),  # Y-axis title style
              axis.title.x = element_text(size = 10, color = "black", angle = 0),   # X-axis title style
              strip.text.x = element_text(size = 10, color = "black"),  # Facet x-label style (not used here)
              strip.text.y = element_text(size = 10, color = "black"),  # Facet y-label style (not used here)
              plot.margin = unit(rep(0.1, 4), 'cm')   # Narrow plot margins (top/right/bottom/left) to save space
            ))

# --- Load and Preprocess Biogeographic & Network Data ---
# Load the main dataset containing network metrics (Connectance, Modularity, Nestedness) and biogeographic predictors (area, isolation)
network_dat <- read_csv("data/network_data.csv")

# Apply log10 transformations to biogeographic predictors to reduce skewness
# This helps meet the linearity assumption of glm (used in stat_smooth later)
network_dat$area <- log10(network_dat$area)  # Log10-transform island area
network_dat$isolation <- log10(network_dat$isolation / 1000)  # Log10-transform isolation (convert to km first)


# --- Create Data Subsets for Comparison (Control vs. Invaded) ---
# Subset 1: Control group (non-invaded sites, invasion == "N"), labeled as "Control"
network_control <- network_dat %>% 
  filter(invasion == "N") %>% 
  mutate(type = "Control")
# Subset 2: Invaded group (all sites reclassified as "Invaded" for comparison with Control)
# Note: `network_inv` is created but not used (likely a leftover from data exploration)
network_inv <- network_dat %>% filter(invasion == "Y")
network_addinv <- network_dat %>% mutate(type = "Invaded")  # All sites labeled as "Invaded"
# Combine Control and Invaded subsets into a single dataframe for group-wise fitting
network_all <- rbind(network_control, network_addinv)
# Set factor levels for "type" to ensure Control appears first in legends/plots
network_all$type <- factor(network_all$type, levels = c("Control", "Invaded"))


# --- Figure 6A: Isolation vs. Connectance ---
# Scatter plot of log10(isolation) vs. Connectance (z-scores), with group-wise glm fits
fig6a <- ggplot() +
  # Scatter points: All sites labeled as "Invaded" (network_addinv), colored by actual invasion status
  geom_point(
    data = network_addinv,
    aes(x = isolation, y = Connectance, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10  # Shape 21 = filled circle with border; alpha = transparency
  ) +
  # Axis labels: Y = Connectance (z-scores); X = no label (shared label in fig6B)
  labs(y = "Connectance (z-scores)", x = "") +
  scale_fill_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  scale_color_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  # Group-wise glm fit lines: Control vs. Invaded (from network_all)
  stat_smooth(
    data = network_all,
    aes(x = isolation, y = Connectance, color = type, fill = type),
    method = "glm", size = 1.5, linetype = 1  # Linetype 1 = solid line
  ) +
  # Add regression equation, R², and p-value (formatted for readability)
  stat_poly_eq(
    data = network_all,
    aes(
      x = isolation, y = Connectance, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~")  # ~~ = line break in plotmath
    ),
    size = 3  # Text size for annotations
  )


# --- Figure 6B: Isolation vs. Modularity ---
# Scatter plot of log10(isolation) vs. Modularity (z-scores), with group-wise glm fits
# Shares x-axis label with fig6A/5C (isolation)
fig6b <- ggplot() +
  geom_point(
    data = network_addinv,
    aes(x = isolation, y = Modularity, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  # Axis labels: Y = Modularity (z-scores); X = Isolation (log10-transformed km)
  labs(y = "Modularity (z-scores)", x = "Distance to closest mainland [log10(km)]") +
  scale_fill_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  scale_color_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  stat_smooth(
    data = network_all,
    aes(x = isolation, y = Modularity, color = type, fill = type),
    method = "glm", size = 1.5, linetype = 1  # Solid line (matches fig6A)
  ) +
  stat_poly_eq(
    data = network_all,
    aes(
      x = isolation, y = Modularity, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~")
    ),
    size = 3
  )


# --- Figure 6C: Isolation vs. Nestedness ---
# Scatter plot of log10(isolation) vs. Nestedness (z-scores), with group-wise glm fits
fig6c <- ggplot() +
  geom_point(
    data = network_addinv,
    aes(x = isolation, y = Nestedness, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Nestedness (z-scores)", x = "")  # X = no label (shared with fig6B)
  scale_fill_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  scale_color_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  stat_smooth(
    data = network_all,
    aes(x = isolation, y = Nestedness, color = type, fill = type),
    method = "glm", size = 1.5, linetype = 1  # Solid line (isolation-related plots)
  ) +
  stat_poly_eq(
    data = network_all,
    aes(
      x = isolation, y = Nestedness, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~")
    ),
    size = 3
  )


# --- Figure 6D: Area vs. Connectance ---
# Scatter plot of log10(area) vs. Connectance (z-scores), with group-wise glm fits
# Uses dashed lines (linetype = 2) to distinguish from isolation-related plots
fig6d <- ggplot() +
  geom_point(
    data = network_addinv,
    aes(x = area, y = Connectance, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Connectance (z-scores)", x = "")  # X = no label (shared with fig6E)
  scale_fill_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  scale_color_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  stat_smooth(
    data = network_all,
    aes(x = area, y = Connectance, color = type, fill = type),
    method = "glm", size = 1.5, linetype = 2  # Linetype 2 = dashed line (area-related plots)
  ) +
  stat_poly_eq(
    data = network_all,
    aes(
      x = area, y = Connectance, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~")
    ),
    size = 3
  )


# --- Figure 6E: Area vs. Modularity ---
# Scatter plot of log10(area) vs. Modularity (z-scores), with group-wise glm fits
# Shares x-axis label with fig6D/5F (area)
fig6e <- ggplot() +
  geom_point(
    data = network_addinv,
    aes(x = area, y = Modularity, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  # Axis labels: Y = Modularity (z-scores); X = Island area (log10-transformed km²)
  labs(y = "Modularity (z-scores)", x = "Island area [log10(km2)]") +
  scale_fill_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  scale_color_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  stat_smooth(
    data = network_all,
    aes(x = area, y = Modularity, color = type, fill = type),
    method = "glm", size = 1.5, linetype = 2  # Dashed line (area-related plots)
  ) +
  stat_poly_eq(
    data = network_all,
    aes(
      x = area, y = Modularity, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~")
    ),
    size = 3
  )


# --- Figure 6F: Area vs. Nestedness ---
# Scatter plot of log10(area) vs. Nestedness (z-scores), with group-wise glm fits
fig6f <- ggplot() +
  geom_point(
    data = network_addinv,
    aes(x = area, y = Nestedness, color = invasion, fill = invasion),
    size = 3, shape = 21, alpha = 5/10
  ) +
  labs(y = "Nestedness (z-scores)", x = "")  # X = no label (shared with fig6E)
  scale_fill_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  scale_color_manual(values = c("#1D91C0", "#FF8000", "#1D91C0", "#FF8000")) +
  stat_smooth(
    data = network_all,
    aes(x = area, y = Nestedness, color = type, fill = type),
    method = "glm", size = 1.5, linetype = 2  # Dashed line (area-related plots)
  ) +
  stat_poly_eq(
    data = network_all,
    aes(
      x = area, y = Nestedness, color = type, fill = type,
      label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~")
    ),
    size = 3
  )


# --- Assemble Final Figure 6 ---
# Combine all 6 subfigures (5A-5F) into a 3-column grid with subfigure labels (A-F)
fig6 <- cowplot::plot_grid(
  fig6a, fig6b, fig6c, fig6d, fig6e, fig6f,
  ncol = 3,  # Arrange plots in 3 columns (rows = 2)
  labels = c('A', 'B', 'C', 'D', 'E', 'F')  # Label each subfigure (top-left corner)
  # Original comment: Extra labels ("G", "H", "I") are commented out (not needed here)
)


# --- Save Final Figure to File ---
# Save as PNG with high resolution (suitable for publications)
# Width = 9 inches (accommodates 3 columns), Height = 6 inches (fits 2 rows)
ggsave(
  filename = paste0('figure/figure6_network_biogeography.png'),  # Output path/name
  plot = fig6,                                               # Plot to save
  width = 9,                                                # Width in inches
  height = 6                                                # Height in inches
)
