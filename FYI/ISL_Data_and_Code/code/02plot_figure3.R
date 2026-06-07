# Load required R packages
library(tidyverse)    # For data manipulation and visualization (includes dplyr, ggplot2, etc.)
library(plot3D)       # For creating 3D plots
library(cowplot)      # For combining multiple ggplot2 plots into a single figure
library(bipartite)    # For bipartite network analysis and visualization
library(ggpubr)       # For creating publication-ready ggplot2 plots and adding statistical tests

# The goal of this script is to plot Figure 3, which compares the predator-prey network structure
# and diet composition between invaded and control ponds.

# Load the network metrics data for all sites/ponds
network_dat <- read_csv("data/network_data.csv")

# Set a global ggplot2 theme for consistent formatting across all plots
# This theme removes grid lines, sets consistent text sizes and colors, and adjusts margins.
theme_set(theme_test() +
            theme(panel.grid = element_blank(),
                  legend.position = "none",          # Remove legend
                  legend.title = element_blank(),
                  legend.text = element_text(color = "black", size = 10),
                  axis.text.x = element_text(size = 10, color = "black"),
                  axis.text.y = element_text(size = 10, color = "black"),
                  axis.title.y = element_text(size = 10, color = "black", angle = 90),
                  axis.title.x = element_text(size = 10, color = "black", angle = 0),
                  strip.text.x = element_text(size = 10, color = "black"), # Facet labels
                  strip.text.y = element_text(size = 10, color = "black"),
                  plot.margin = unit(rep(0.1, 4), 'cm'))) # Reduce plot margins

# --- Plot Figure 3a: 3D Scatter Plot of Network Metrics ---

# Prepare data for the 3D plot, filtering out rows where Connectance is NA
plot3d_netw <- network_dat %>% filter(!is.na(Connectance))

# Convert the 'invasion' column to a factor for proper color mapping
plot3d_netw$invasion <- as.factor(plot3d_netw$invasion)

# Open a CairoPNG device to save the 3D plot with high resolution
# This is used instead of png() for better anti-aliasing and font rendering.
Cairo::CairoPNG(
  filename = paste0("figure/figure3a_network_3d.png "), # Output file path
  width = 10, height = 10, units = "in",              # Size in inches
  dpi = 600)                                          # High resolution for publication

# Create a 3D scatter plot using plot3D::scatter3D
scatter3D(
  x = plot3d_netw$Connectance,   # X-axis: Connectance values
  y = plot3d_netw$Nestedness,    # Y-axis: Nestedness values
  z = plot3d_netw$Modularity,    # Z-axis: Modularity values
  colvar = as.numeric(plot3d_netw$invasion), # Color points by invasion status (control vs. invaded)
  col = c("#1D91C0", "#FF8000"),             # Colors for the two groups (blue for control, orange for invaded)
  pch = 16,                                  # Point type (solid circle)
  theta = 30,                                # Azimuthal angle (rotation around z-axis)
  phi = 20,                                  # Polar angle (elevation)
  bty = "b2",                                # Type of box around the plot
  xlab = "Connectance",                      # X-axis label
  ylab = "Nestedness",                       # Y-axis label
  zlab = "Modularity",                       # Z-axis label
  main = "Network Structure"                 # Plot title
)

# Close the PNG device and save the file
dev.off()


# --- Plot Figure 3a Insets: Example Network Webs ---

# --- Example Site 1 (control pond) ---

# Load the interaction data for the first example site
ex_site1 <- read_csv("data/example_site1_JT_site060.csv")

# Reshape the data into a matrix (prey as rows, predators as columns)
# using RRA (Relative Read Abundance) as the interaction strength.
ex_site1_matrix <- ex_site1 %>%
  reshape2::acast(formula = prey ~ species, fun.aggregate = sum, value.var = "RRA") %>%
  as.data.frame()

# Define the sequence (order) of species for plotting the network web.
# This is done to arrange the species in a meaningful order (e.g., by dominance).
seq.high <- colnames(ex_site1_matrix) # Order of predator species (top level)

# Determine the order of prey species (bottom level) based on their maximum RRA with any predator.
# This helps in arranging the prey to visually highlight strong interactions.
site_matrix_prec <- ex_site1 %>%
  group_by(species) %>%
  dplyr::arrange(desc(RRA)) %>%       # Arrange within each predator from highest to lowest RRA
  ungroup() %>%
  dplyr::arrange(match(species, seq.high)) %>% # Reorder rows to match the predator sequence
  group_by(prey) %>%
  filter(RRA == max(RRA))             # For each prey, keep only its strongest interaction

seq.low <- unique(site_matrix_prec$prey) # Order of prey species (bottom level)
seq_spec <- list(seq.high = seq.high, seq.low = seq.low) # Combine into a list for plotweb

# Define colors for the network web (all blue for this control site example)
high.colors <- c("#1D91C0")
low.colors <- c("#1D91C0")
int.colors <- c("#1D91C0")

# Open a CairoPNG device to save the first example network web
Cairo::CairoPNG(
  filename = paste0("figure/figure3a_ex1.png "),
  width = 8, height = 5, units = "in", dpi = 600)

# Plot the bipartite network web using bipartite::plotweb
plotweb(ex_site1_matrix,
        method = "normal",              # Plotting method
        arrow = "center",               # Arrow style
        col.low = low.colors,           # Color for prey (bottom) nodes
        col.high = high.colors,         # Color for predator (top) nodes
        col.interaction = (int.colors), # Color for interaction arrows
        sequence = seq_spec,            # Order of species on both levels
        bor.col.low = "transparent",    # Transparent border for prey nodes
        bor.col.high = "transparent",   # Transparent border for predator nodes
        bor.col.interaction = "transparent") # Transparent border for arrows

# Close the PNG device
dev.off()


# --- Example Site 2 (invaded pond) ---

# Load the interaction data for the second example site
ex_site2 <- read_csv("data/example_site2_DS_site021.csv")

# Reshape the data into a matrix
ex_site2_matrix <- ex_site2 %>%
  reshape2::acast(formula = prey ~ species, fun.aggregate = sum, value.var = "RRA") %>%
  as.data.frame()

# Define the sequence of species for plotting (same method as for site 1)
seq.high <- colnames(ex_site2_matrix)
site_matrix_prec <- ex_site2 %>%
  group_by(species) %>%
  dplyr::arrange(desc(RRA)) %>%
  ungroup() %>%
  dplyr::arrange(match(species, seq.high)) %>%
  group_by(prey) %>%
  filter(RRA == max(RRA))
seq.low <- unique(site_matrix_prec$prey)
seq_spec <- list(seq.high = seq.high, seq.low = seq.low)

# Define colors for the network web (all orange for this invaded site example)
high.colors <- c("#FF8000")
low.colors <- c("#FF8000")
int.colors <- c("#FF8000")

# Open a CairoPNG device to save the second example network web
Cairo::CairoPNG(
  filename = paste0("figure/figure3a_ex2.png "),
  width = 8, height = 5, units = "in", dpi = 600)

# Plot the bipartite network web
plotweb(ex_site2_matrix,
        method = "normal", arrow = "center",
        col.low = low.colors, col.high = high.colors, col.interaction = (int.colors),
        sequence = seq_spec,
        bor.col.low = "transparent", bor.col.high = "transparent", bor.col.interaction = "transparent")

# Close the PNG device
dev.off()


# --- Plot Figure 3b-d: Boxplots of Network Metrics ---

# Create a boxplot for Connectance, grouped by invasion status
fig3b <- network_dat %>%
  ggplot(aes(x = invasion, y = Connectance, fill = invasion)) +
  geom_boxplot(width = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c("#3498DB", "#FF8000")) + # Colors for control (blue) and invaded (orange)
  scale_color_manual(values = c("#3498DB", "#FF8000")) +
  xlab("") + ylab("Connectance (z-scores)") + # Set axis labels
  # Add Kruskal-Wallis test result as text on the plot
  stat_kruskal_test(label = "{p.format}{p.signif}", p.adjust.method = "holm")

# Create a boxplot for Modularity, grouped by invasion status
fig3c <- network_dat %>%
  ggplot(aes(x = invasion, y = Modularity, fill = invasion)) +
  geom_boxplot(width = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c("#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000")) +
  xlab("") + ylab("Modularity (z-scores)") +
  stat_kruskal_test(label = "{p.format}{p.signif}", p.adjust.method = "holm")

# Create a boxplot for Nestedness, grouped by invasion status
fig3d <- network_dat %>%
  ggplot(aes(x = invasion, y = Nestedness, fill = invasion)) +
  geom_boxplot(width = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c("#3498DB", "#FF8000")) +
  scale_color_manual(values = c("#3498DB", "#FF8000")) +
  xlab("") + ylab("Nestedness (z-scores)") + # Corrected the typo here
  stat_kruskal_test(label = "{p.format}{p.signif}", p.adjust.method = "holm")

# Combine the three boxplots into a single figure using cowplot::plot_grid
fig3df <- plot_grid(fig3b, fig3c, fig3d,
                    ncol = 3,             # Arrange plots in 3 columns
                    labels = c('B', 'C', 'D')) # Add labels to each subplot

# Save the combined figure to a file
ggsave(filename = paste0('figure/figure3b-d_network_change.png'),
       plot = fig3df,
       width = 9, height = 3)
```
