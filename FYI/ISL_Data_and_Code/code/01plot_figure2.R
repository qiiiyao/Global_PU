```r
# Load required R packages for data manipulation, scientific visualization, and plot styling
library(tidyverse)    # Core package for data cleaning (dplyr) and visualization (ggplot2)
library(ggsci)        # For scientific journal-style color palettes
library(paletteer)    # For accessing pre-defined color palettes (used for prey/species coloring)
library(cartography)  # For spatial data visualization (not explicitly used here, likely a leftover from prior analyses)
library(ggpubr)       # For creating publication-ready plots (sets clean themes and formatting)

# Read dataset for Figure 2a (prey composition bar plot)
plot_dat <- read_csv("data/figs2_data.csv")

# Specify the order of prey taxa (x-axis grouping in bar plot) as a factor
plot_dat$prey <- factor(plot_dat$prey, levels = c("Isopoda", "Hymenoptera", "Diptera", "Hemiptera", "Coleoptera", 
                                                  "Lepidoptera", "Araneae", "Sarcoptiformes", "Thysanoptera", 
                                                  "Trombidiformes", "Dermaptera", "Polydesmida", "Orthoptera", 
                                                  "Podocopida", "Decapoda", "Pseudoscorpiones", "Poduromorpha", 
                                                  "Entomobryomorpha", "Mesostigmata", "Diplostraca", "Odonata", 
                                                  "Trichoptera", "Blattodea", "Amphipoda", "Scolopendromorpha", 
                                                  "Phthiraptera", "Halocyprida", "Geophilomorpha", "Symphypleona", 
                                                  "Calanoida", "Cyclopoida", "Ephemeroptera", "Psocoptera"))

# Specify the order of amphibian species (x-axis in bar plot) as a factor
plot_dat$species <- factor(plot_dat$species, levels = c("Fejervarya multistriata", "Microhyla fissipes",
                                                        "Pelophylax nigromaculatus", "Bufo gargarizans",
                                                        "Pelophylax plancyi", "Hyla chinensis", 
                                                        "Rana zhenhaiensis", "Polypedates braueri",
                                                        "Quasipaa spinosa", "Hylarana latouchii",  
                                                        "Aquarana catesbeianus"))  # Invasive American bullfrog (last for emphasis)


# --- Figure 2a: Stacked Bar Plot of Prey Composition (RRA) by Species ---
# Visualizes the relative contribution of each prey taxon to the diet of each amphibian species
figs2a <- ggplot(plot_dat, aes(x = species, y = RRA)) +
  # Stacked bar plot: each bar represents a species, segments represent prey taxa (RRA sum to 1)
  geom_bar(aes(color = prey, fill = prey), stat = "identity") +
  # # Commented-out code: Likely for adding count labels to bar segments (not used here)
  # geom_text(aes(y = lab_pos, label = counts, group = color), size = 3) +
  # Set publication-style theme: legend on the right, x-axis labels rotated 45° to avoid overlap
  theme_pubr(legend = "right", x.text.angle = 45) +
  # Use discrete color palette from ggsci ("default_igv") for prey taxa (consistent with scientific standards)
  scale_color_paletteer_d("ggsci::default_igv", direction = 1) +
  scale_fill_paletteer_d("ggsci::default_igv", direction = 1)  # Match fill color to border color
  # # Commented-out code: Likely for faceting by an additional variable (e.g., "fly")
  # + facet_wrap(~fly)


# --- Prepare Data for Figure 2b (NMDS Ordination Plot) ---
# Read results of Non-metric Multidimensional Scaling (NMDS) analysis
# NMDS summarizes prey composition dissimilarity between species in 2-dimensional space
nmds_result <- read_csv("data/figs2_nmds_result.csv")

# Specify the order of species for NMDS plot (consistent with Figure 2a)
nmds_result$species <- factor(nmds_result$species, levels = c("Fejervarya_multistriata", 
                                                              "Pelophylax_nigromaculatus", "Microhyla_fissipes", 
                                                              "Bufo_gargarizans", "Rana_zhenhaiensis", "Pelophylax_plancyi",
                                                              "Hylarana_latouchii", "Polypedates_braueri", 
                                                              "Hyla_chinensis", "Aquarana_catesbeianus"))  # Invasive bullfrog last

# Define species-specific point shapes (ensures each species has a unique symbol)
shape_list <- c(
  "Fejervarya_multistriata" = 16,    # Solid circle
  "Microhyla_fissipes" = 0,          # Open square
  "Pelophylax_nigromaculatus" = 1,   # Open circle
  "Pelophylax_plancyi" = 8,          # Octagon
  "Bufo_gargarizans" = 3,            # Plus sign
  "Rana_zhenhaiensis" = 4,           # X sign
  "Hylarana_latouchii" = 5,          # Diamond
  "Hyla_chinensis" = 6,              # Triangle (point up)
  "Polypedates_megacephalus" = 7,    # Square with cross
  "Aquarana_catesbeianus" = 17     # Solid triangle (point up, for invasive bullfrog)
)

# Define species-specific colors (ensures each species has a unique color; invasive bullfrog = red)
color_list <- list(
  "Fejervarya_multistriata" = "#41B7C4",    # Teal
  "Microhyla_fissipes" = "#5CC0C0",          # Light teal
  "Pelophylax_nigromaculatus" = "#7ECDBB",   # Mint green
  "Bufo_gargarizans" = "#98D7B7",           # Pale green
  "Rana_zhenhaiensis" = "#B6E4B3",           # Very pale green
  "Pelophylax_plancyi" = "#C7EBB1",          # Almost yellow-green
  "Hylarana_latouchii" = "#EDF8BC",          # Pale yellow
  "Hyla_chinensis" = "#F1FABF",              # Light yellow
  "Polypedates_megacephalus" = "#FCFED4",    # Near white
  "Aquarana_catesbeianus" = "#DC0000"      # Bright red (highlights invasive species)
)

# NMDS stress value: Measures how well the 2D ordination represents the original dissimilarity matrix
# Stress = 0.149: Acceptable (stress < 0.2 indicates a useful ordination)
stress <- "0.149"


# --- Figure 2b: NMDS Ordination Plot of Species Prey Composition ---
# Visualizes dissimilarity in prey composition between amphibian species
figs2b <- ggplot(data = nmds_result, aes(x = MDS1, y = MDS2,
                                         color = species, fill = species,
                                         shape = species)) +
  # Scatter points: Each point represents a species; size = 4 for visibility
  geom_point(size = 4) +
  # # Commented-out code: Adds convex hulls around species (not used; replaced with confidence ellipses)
  # stat_chull(alpha = 0.1, geom = "polygon") +
  # Axis labels and caption (stress value for transparency)
  labs(x = paste("NMDS 1"),  # X-axis: First NMDS dimension
       y = paste("NMDS 2"),  # Y-axis: Second NMDS dimension
       caption = paste('Stress =', stress))  # Caption: NMDS stress value
  # Set custom theme for clarity and professionalism
  theme(
    # legend.position = c(0.9, 0.8),  # Commented-out: Alternative legend position (top-right)
    legend.title = element_blank(),          # Remove legend title (species are self-explanatory)
    panel.grid = element_blank(),            # Remove all grid lines (reduces clutter)
    plot.title = element_text(hjust = 0),    # Left-align plot title (if used)
    panel.grid.major = element_blank(),      # Confirm no major grid lines
    panel.grid.minor = element_blank(),      # Confirm no minor grid lines
    panel.background = element_rect(color = 'black', fill = 'transparent'),  # Black border, transparent background
    axis.text = element_text(color = "black", size = 10)  # Black axis text (size = 10 for readability)
  ) +
  # Add dashed gray lines at x=0 and y=0 (reference lines for ordination origin)
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")


# Add 90% confidence ellipses to the NMDS plot
# Ellipses represent the expected range of variation for each species' prey composition
figs2b <- figs2b + stat_ellipse(
  data = nmds_result,
  geom = "polygon",          # Draw ellipses as filled polygons
  level = 0.9,               # 90% confidence level
  linetype = 2,              # Dashed ellipse borders
  linewidth = 0.5,           # Thin border lines
  aes(fill = species),       # Fill ellipses with species-specific colors
  alpha = 0.3,               # Transparency (alpha = 0.3) to avoid obscuring points
  show.legend = T            # Show ellipses in the legend (matches points)
) +
  # Re-apply species-specific fill colors (ensures consistency with points)
  scale_fill_manual(values = color_list)


# --- Assemble and Save Combined Figure (Figure S1) ---
# Combine Figure 2a (bar plot) and Figure 2b (NMDS plot) into a single vertical stack (1 column)
# Add labels "A" and "B" to identify subfigures (standard for scientific publications)
figs2 <- cowplot::plot_grid(
  figs2a, figs2b,
  ncol = 1,                  # Arrange plots vertically (1 column, 2 rows)
  labels = c('A', 'B')       # Label subfigures (top-left corner of each plot)
)

# --- Figure 2c: Vertebrate prey composition in the diet of the invasive bullfrog ---
bullfrog_vert_diet <- read_csv("data/fig2_bullfrog_vert_diet.csv")
bvdt <- bullfrog_vert_diet %>%
  arrange(desc(RRA)) %>%
  mutate(
    percent = round(RRA * 100, 1),
    label = ifelse(percent >= 5, paste0(percent, "%"), ""),
    ypos = cumsum(RRA) - 0.5 * RRA
  )
bvdt$Prey <- factor(bvdt$Prey,levels =c(bvdt$Prey))
fig2c <- ggplot(bvdt, aes(x = "", y = RRA, fill = Prey)) +
  geom_bar(stat = "identity", width = 0.3, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = ypos, label = label), 
            color = "black", size = 3) +
  scale_fill_brewer(palette = "Set3") +
  # theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )
fig2c
fig2bc <- cowplot::plot_grid(
  fig2b,fig2c,
  ncol = 2, labels = c('B', 'C'))  
fig2 <- cowplot::plot_grid(
  fig2a,fig2bc,
  ncol = 1, labels = c('A')) 
fig2
ggsave(paste0("figure/figure2.pdf"), fig2, width = 6.5, height = 6.5)


