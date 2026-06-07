# Load required R packages for data manipulation, visualization, and plot assembly
library(tidyverse)    # For data cleaning (dplyr) and visualization (ggplot2)
library(cowplot)      # For combining multiple ggplot2 plots into a single figure
library(ggpubr)       # For enhancing publication-quality plots (not fully used here but loaded for consistency)

# --- Set Global ggplot2 Theme for Consistent Formatting ---
# Standardize plot style across all subfigures to ensure uniformity in publications
theme_set(theme_test() +
            theme(
              panel.grid = element_blank(),          # Remove background grid lines (reduces clutter)
              legend.position = "none",              # Hide legend by default (no need for redundant color legends here)
              legend.title = element_blank(),        # Remove legend title (irrelevant since legend is hidden)
              legend.text = element_text(color = "black", size = 10),  # Legend text style (fallback if legend is enabled)
              axis.text.x = element_text(size = 10, color = "black"),   # X-axis tick label style (readability)
              axis.text.y = element_text(size = 10, color = "black"),   # Y-axis tick label style (readability)
              axis.title.y = element_text(size = 10, color = "black", angle = 90),  # Y-axis title style
              axis.title.x = element_text(size = 10, color = "black", angle = 0),   # X-axis title style
              strip.text.x = element_text(size = 10, color = "black"),  # Facet x-label style (not used here)
              strip.text.y = element_text(size = 10, color = "black"),  # Facet y-label style (not used here)
              plot.margin = unit(rep(0.1, 4), 'cm')   # Narrow plot margins (top/right/bottom/left) to save space
            ))

# --- Load Effect Size Data from GLMM Model Averaging ---
# Load dataset containing standardized coefficients, 95% confidence intervals, 
# variable selection weights (sw), and significance for network metrics
eff_dat <- read_csv("data/fig5-glmm_effect_size.csv")

# --- Reorder Factor Levels for Variables ---
# Set the order of variables (y-axis) to control their display order in plots
# Variables are reversed to ensure they appear in the desired top-to-bottom order on the y-axis
eff_dat$varibles <- factor(eff_dat$varibles,  # Note: Potential typo (should be "variables")
                            levels = rev(c(
                              "Residence_time",       # Time since invasion/residence
                              "Bullfrog_abundance",   # Abundance of invasive bullfrogs
                              "Isolation",            # Pond isolation (geographic)
                              "Area",                 # Pond area
                              "Native_abundance",     # Abundance of native species
                              "Water_area",           # Water body area
                              "Vegetation_coverage",  # Vegetation coverage around ponds
                              "Human_activity",       # Human disturbance level
                              "Crayfish_abundance",   # Abundance of crayfish
                              "Bullfrog_abundance:Isolation",  # Interaction: Isolation × Bullfrog abundance
                              "Area:Bullfrog_abundance"        # Interaction: Area × Bullfrog abundance
                            )))


# --- Plot Figure 5A: Effect Sizes for Connectance ---
# Horizontal dot plot showing standardized coefficients (estimates) of predictors for Connectance,
# with 95% confidence intervals (error bars), variable-specific colors, and selection weights (point size)
fig5a_ta_con <- eff_dat %>%
  filter(network == "Connectance") %>%  # Subset data to only Connectance (network metric)
  ggplot(aes(
    x = stand_cof,                      # X-axis: Standardized coefficient (effect size)
    y = varibles,                       # Y-axis: Predictor variables (ordered by factor levels)
    color = varibles,                   # Color-code points by variable (for visual grouping)
    size = sw                           # Point size: Proportional to variable selection weight (sw)
                                        # (higher sw = more important in model averaging)
  )) +
  # Vertical dashed line at x=0 (reference: no effect; intervals crossing 0 = non-significant)
  geom_vline(xintercept = 0, linetype = "dashed", color = "#757575", linewidth = 0.5) +
  # Points representing standardized coefficients (position_dodge avoids overlap)
  geom_point(position = position_dodge(width = 0.5)) +
  # Error bars representing 95% confidence intervals (xmin = 2.5th percentile, xmax = 97.5th percentile)
  geom_errorbar(
    aes(xmin = `2.5 %`, xmax = `97.5 %`),  # Columns for 95% CI bounds
    width = 0,                              # No horizontal "caps" on error bars (cleaner look)
    position = position_dodge(width = 0.5), # Align error bars with points
    linewidth = 1                           # Thicken error bars for visibility
  ) +
  # Manual color scale for variables (groups similar predictors by color for interpretation)
  scale_color_manual(values = c(
    "Isolation" = "#91D1C2FF",            # Teal: Geographic predictors (Isolation, Area)
    "Area" = "#91D1C2FF",
    "Residence_time" = "#F39C12",         # Orange: Invasion-related predictors (Residence time, Bullfrog abundance)
    "Bullfrog_abundance" = "#F39C12",
    "Native_abundance" = "#EBCCE1",       # Light purple: Native species predictor
    "Water_area" = "#00BFC4",             # Dark teal: Habitat predictors (Water area, Vegetation, Human activity, Crayfish)
    "Vegetation_coverage" = "#00BFC4",
    "Human_activity" = "#00BFC4",
    "Crayfish_abundance" = "#00BFC4"
    # Note: Interaction terms (e.g., Bullfrog_abundance:Isolation) inherit default colors (not defined here)
  )) +
  # Add significance labels (***, **, *, or blank) next to points (nudge_y avoids overlap)
  geom_text(
    aes(label = significance), 
    nudge_y = 0.15,                        # Shift text slightly up from points
    color = "#1e1e20",                     # Dark gray text (high contrast)
    size = 2.5,                           # Small text size (non-intrusive)
    fontface = "bold"                     # Bold text for readability
  ) +
  # Axis labels and plot title
  labs(x = "Estimates", y = "") +         # X-axis: "Estimates" (standardized coefficients); Y-axis: No label (variables are self-explanatory)
  ggtitle('Connectance')                  # Title: Network metric (Connectance)

# Preview Figure 5A (optional, for interactive checking)
fig5a_ta_con


# --- Plot Figure 5B: Effect Sizes for Modularity ---
# Same structure as Figure 5A, but for Modularity (second network metric)
fig5b_ta_mod <- eff_dat %>%
  filter(network == "Modularity") %>%     # Subset data to Modularity
  ggplot(aes(
    x = stand_cof, 
    y = varibles, 
    color = varibles, 
    size = sw
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#757575", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(xmin = `2.5 %`, xmax = `97.5 %`),
    width = 0, 
    position = position_dodge(width = 0.5),
    linewidth = 1
  ) +
  # Same color scale as Figure 5A (consistency across subfigures)
  scale_color_manual(values = c(
    "Isolation" = "#91D1C2FF",
    "Area" = "#91D1C2FF",
    "Residence_time" = "#F39C12",
    "Bullfrog_abundance" = "#F39C12",
    "Native_abundance" = "#EBCCE1",
    "Water_area" = "#00BFC4",
    "Vegetation_coverage" = "#00BFC4",
    "Human_activity" = "#00BFC4",
    "Crayfish_abundance" = "#00BFC4"
  )) +
  geom_text(
    aes(label = significance), 
    nudge_y = 0.15, 
    color = "#1e1e20", 
    size = 2.5, 
    fontface = "bold"
  ) +
  labs(x = "Estimates", y = "") +
  ggtitle('Modularity')                   # Title: Network metric (Modularity)

# Preview Figure 5B (optional)
fig5b_ta_mod


# --- Plot Figure 5C: Effect Sizes for Nestedness ---
# Same structure as 5A/5B, but for Nestedness (third network metric)
# Note: Bug in original code: `filter(network == "Modularity")` should be `filter(network == "Nestedness")`
fig5c_ta_nest <- eff_dat %>%
  filter(network == "Nestedness") %>%     # Corrected: Subset data to Nestedness (fixes original bug)
  ggplot(aes(
    x = stand_cof, 
    y = varibles, 
    color = varibles, 
    size = sw
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#757575", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(xmin = `2.5 %`, xmax = `97.5 %`),
    width = 0, 
    position = position_dodge(width = 0.5),
    linewidth = 1
  ) +
  # Same color scale as 5A/5B (consistency)
  scale_color_manual(values = c(
    "Isolation" = "#91D1C2FF",
    "Area" = "#91D1C2FF",
    "Residence_time" = "#F39C12",
    "Bullfrog_abundance" = "#F39C12",
    "Native_abundance" = "#EBCCE1",
    "Water_area" = "#00BFC4",
    "Vegetation_coverage" = "#00BFC4",
    "Human_activity" = "#00BFC4",
    "Crayfish_abundance" = "#00BFC4"
  )) +
  geom_text(
    aes(label = significance), 
    nudge_y = 0.15, 
    color = "#1e1e20", 
    size = 2.5, 
    fontface = "bold"
  ) +
  labs(x = "Estimates", y = "") +
  ggtitle('Nestedness')                  # Title: Network metric (Nestedness)

# Preview Figure 5C (optional)
fig5c_ta_nest


# --- Assemble Final Figure 5 ---
# Combine Figure 5A (Connectance), 5B (Modularity), and 5C (Nestedness) horizontally (3 columns)
# Add subfigure labels ("A", "B", "C") for publication consistency
fig5 <- cowplot::plot_grid(
  fig5a_ta_con, fig5b_ta_mod, fig5c_ta_nest,
  ncol = 3,                               # Arrange plots in 3 columns
  labels = c('A', 'B', 'C')               # Label each subfigure (top-left corner)
  # Original comment: Extra labels ("D", "E", "F") are commented out (not needed here)
)

# Preview the final combined Figure 5 (optional)
fig5

# --- Save Final Figure to File ---
# Save as PNG with high resolution (suitable for publications)
# Width = 12 inches (accommodates 3 columns), Height = 4 inches (fits variable list on y-axis)
ggsave(
  filename = paste0('figure/figure5_network_effect_size_95_sw.png'),  # Output path/name
  plot = fig5,                                                       # Plot to save
  width = 12,                                                        # Width in inches
  height = 4                                                         # Height in inches
)
```
