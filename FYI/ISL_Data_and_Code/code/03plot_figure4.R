# Load required R packages for data manipulation, visualization, and plot assembly
library(tidyverse)    # For data cleaning, manipulation (dplyr) and visualization (ggplot2)
library(ggpubr)       # For creating publication-quality plots and adding statistical annotations
library(ggpmisc)      # For enhancing ggplot2 with additional statistical elements (e.g., test labels)
library(cowplot)      # For combining multiple ggplot2 plots into a single cohesive figure

# --- Load Data for Centrality Analysis ---
# Load dataset containing species-level centrality metrics (degree and closeness) 
# and invasion status for each site/pond
centrality_dat <- read_csv("data/fig4_centrality_data.csv")

# --- Set Global ggplot2 Theme for Consistent Formatting ---
# Standardize plot style across all figures: remove grid lines, set legend/font properties, adjust margins
theme_set(theme_test() +
            theme(panel.grid = element_blank(),          # Remove background grid lines
                  legend.position = "none",              # Hide legend by default (adjust per plot if needed)
                  legend.title = element_blank(),        # Remove legend title
                  legend.text = element_text(color = "black", size = 10),  # Legend text style
                  axis.text.x = element_text(size = 10, color = "black"),   # X-axis tick label style
                  axis.text.y = element_text(size = 10, color = "black"),   # Y-axis tick label style
                  axis.title.y = element_text(size = 10, color = "black", angle = 90),  # Y-axis title style
                  axis.title.x = element_text(size = 10, color = "black", angle = 0),   # X-axis title style
                  strip.text.x = element_text(size = 10, color = "black"),  # Facet x-label style
                  strip.text.y = element_text(size = 10, color = "black"),  # Facet y-label style
                  plot.margin = unit(rep(0.1, 4), 'cm'))) # Narrow plot margins (top/right/bottom/left)

library(hrbrthemes)   # Optional package for additional theme options (not fully used here, but loaded for potential styling)

# --- Define Point Shapes for Each Amphibian Species ---
# Assign unique plot symbols to distinguish different amphibian species in the scatter plot
# Shape codes follow R's standard plotting symbols (e.g., 16 = solid circle, 0 = open square)
shape_list <- c(
  "Fejervarya_multistriata" = 16,    # Solid circle (Paddy frog)
  "Microhyla_fissipes" = 0,          # Open square (Ornate chorus frog)
  "Pelophylax_nigromaculatus" = 1,   # Open circle (Black-spotted pond frog)
  "Pelophylax_plancyi" = 8,          # Octagon (Gold spotted pond frog)
  "Bufo_gargarizans" = 3,            # Plus sign (Asian toad)
  "Rana_zhenhaiensis" = 4,           # X sign (Zhenhai brown frog)
  "Hylarana_latouchii" = 5,          # Diamond (Brown wood frog)
  "Hyla_chinensis" = 6,              # Triangle (Chinese tree frog)
  "Polypedates_braueri" = 7,    # Square with cross (Brauer's whipping frog)
  "Aquarana_catesbeianus" = 17     # Solid triangle (American bullfrog, invasive species)
)

# --- Plot figure 4a: Scatter Plot of Centrality Metrics ---
# Scatter plot comparing two centrality metrics (Degree vs. Closeness) across species,
# colored by invasion status (control vs. invaded) and shaped by species identity
fig4a <- ggplot(centrality_dat, aes(
  x = degree,                          # X-axis: Degree Centrality (z-scores, standardized)
  y = closeness,                       # Y-axis: Closeness Centrality (z-scores, standardized)
  shape = species,                     # Differentiate species by point shape
  color = invasion                     # Differentiate invasion status by color
)) +  
  geom_point(size = 3) +               # Plot points with size 3 (visible but not overwhelming)
  # Set colors for invasion status: blue = control, orange = invaded
  scale_color_manual(values = c("#3498DB", "#F39C12")) +
  scale_fill_manual(values = c("#3498DB", "#F39C12")) +  # Consistent fill color (redundant here but for consistency)
  # Assign species-specific shapes using the pre-defined `shape_list`
  scale_shape_manual(values = shape_list) +
  # Axis labels (standardized centrality metrics in z-scores)
  xlab("Degree Centrality (z-scores)") + 
  ylab("Closeness Centrality (z-scores)") +
  theme_bw() +                         # Apply black-and-white theme (complements earlier global theme)
  theme(legend.position = "bottom")    # Move legend to bottom to avoid overlapping with data points

# Add subplot label "A" to figure 4a using cowplot (consistent with publication labeling)
fig4a <- cowplot::plot_grid(fig4a, ncol = 1, labels = c('A'))  


# --- Plot figure 4b & 3c: Diet Change Analysis (Boxplots) ---
# Load dataset containing diet metrics (prey number, diet breadth) for native species and bullfrogs
diet_change_dat <- read_csv("data/fig4b-c_diet_change_data.csv")

# Subset data: Isolate data for native amphibian species
# Include only key native species (predefined list) for comparison
diet_change_nat <- diet_change_dat %>% filter(species %in% c(
  "Fejervarya_multistriata", 
  "Pelophylax_nigromaculatus", 
  "Microhyla_fissipes", 
  "Bufo_gargarizans"
))

# Reorder factor levels for native species to control their display order on the x-axis
# Ensures consistent and logical species ordering in figure 4b and 3c
diet_change_nat$species <- factor(diet_change_nat$species, levels = c(
  "Fejervarya_multistriata", 
  "Pelophylax_nigromaculatus", 
  "Microhyla_fissipes", 
  "Bufo_gargarizans"
))


# --- figure 4b: Boxplot of Native Species' Prey Number (log10-transformed) ---
# Compare log10(prey_no) (species degree, i.e., number of prey) of native species 
# between control and invaded ponds
fig4b <- diet_change_nat %>%
  ggplot(aes(
    x = species,                       # X-axis: Native amphibian species
    y = log10(prey_no),                # Y-axis: log10-transformed number of prey (reduces skewness)
    fill = invasion                     # Fill boxplots by invasion status (control vs. invaded)
  )) +
  geom_boxplot(width = 0.5, alpha = 0.8) +  # Boxplots with moderate width and transparency (0.8 = slight overlap allowed)
  # Consistent color scheme: blue = control, orange = invaded
  scale_fill_manual(values = c("#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000")) +
  # Perform Kruskal-Wallis test (non-parametric ANOVA) to compare groups; show p-value and significance stars
  stat_kruskal_test(
    label = "{p.format} {p.signif}",  # Label format: p-value + significance symbols (e.g., *, **)
    p.adjust.method = "holm"          # Holm-Bonferroni correction for multiple comparisons
  ) +
  xlab("") +                          # Remove x-axis title (species names are self-explanatory)
  ylab("The number of prey")              # Y-axis title (matches title for clarity)


# --- figure 4c: Boxplot of Native Species' Standardized Diet Breadth ---
# Compare standardized diet breadth of native species between control and invaded ponds;
# add bullfrog average as reference
fig4c <- diet_change_nat %>%
  ggplot(aes(
    x = species,                       # X-axis: Native amphibian species
    y = scaled.breadth,                # Y-axis: Standardized diet breadth (normalized for comparison)
    fill = invasion                     # Fill boxplots by invasion status
  )) +
  geom_boxplot(width = 0.5, alpha = 0.8) +  # Consistent boxplot style with figure 4b
  scale_fill_manual(values = c("#3498DB", "#FF8000")) +  # Match color scheme
  scale_color_manual(values = c("#3498DB", "#FF8000")) +
  # Add red dashed line: Average standardized diet breadth of invasive bullfrogs (reference)
  # Use `na.rm = TRUE` to ignore any NA values in bullfrog data (avoids errors)
  geom_hline(
    yintercept = mean(bullfrog_m$scaled.breadth, na.rm = TRUE), 
    linetype = "dashed", 
    color = "red", 
    size = 1
  ) +
  # Kruskal-Wallis test with Holm correction (same as figure 4b)
  stat_kruskal_test(
    label = "{p.format}{p.signif}", 
    p.adjust.method = "holm"
  ) +
  xlab("") +                          # Remove x-axis title
  ylab("Diet breadth")                # Y-axis title

# --- Assemble figure 4: Combine 3a, 3b, and 3c ---
# First, stack figure 4b and 3c vertically (1 column) with labels "B" and "C"
fig4bc <- cowplot::plot_grid(fig4b, fig4c, ncol = 1, labels = c('B', 'C'))  

# Then, arrange figure 4a (left) and fig4bc (right) horizontally (2 columns)
fig4 <- cowplot::plot_grid(fig4a, fig4bc, ncol = 2)  

# Save the final combined figure 4 to a file (publication-ready resolution)
# Width = 8.5 inches, Height = 5 inches (balances readability and space)
ggsave(
  filename = paste0('figure/figure4.png'), 
  plot = fig4, 
  width = 8.5, 
  height = 5
)




