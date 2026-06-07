# Load required R packages for data manipulation, visualization, and plot assembly
library(tidyverse)    # For data cleaning (dplyr) and visualization (ggplot2)
library(ggpubr)       # For enhancing publication-quality plots
library(ggpmisc)      # For adding regression equations and statistics to plots
library(cowplot)      # For combining multiple ggplot2 plots into a single figure

# --- Figure S3a: Diet Breadth vs. Isolation (Non-invaded sites only) ---

# Load diet change data and combine it with network data
# Filter to include only non-invaded sites (invasion == "N")
diet_change_dat <- read_csv("data/fig4b-c_diet_change_data.csv") %>%
  left_join(network_dat) %>%  # Assuming 'network_dat' is already loaded in the environment
  filter(invasion == "N")

# Set a global ggplot2 theme for consistent, publication-ready formatting
theme_set(theme_test() +
            theme(
              panel.grid = element_blank(),          # Remove background grid lines
              legend.position = "none",              # Hide legend (not needed for these plots)
              legend.title = element_blank(),        # Remove legend title
              legend.text = element_text(color = "black", size = 10),
              axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.title.y = element_text(size = 10, color = "black", angle = 90),
              axis.title.x = element_text(size = 10, color = "black", angle = 0),
              strip.text.x = element_text(size = 10, color = "black"),
              strip.text.y = element_text(size = 10, color = "black"),
              plot.margin = unit(rep(0.1, 4), 'cm')   # Narrow plot margins
            ))

# Create a scatter plot with a linear regression line for non-invaded sites
figs2a <- ggplot(data = diet_change_dat, aes(x = isolation, y = scaled.breadth)) +
  # Scatter points: all black, semi-transparent
  geom_point(size = 3, shape = 21, alpha = 5/10, fill = "black", colour = "black") +
  # Axis labels
  labs(y = "Diet breadth", x = "Distance to closest mainland [log10(km)]") +
  # Add a linear regression line (method = "lm")
  stat_smooth(method = "lm", size = 1.5, linetype = 1, colour = "black") +
  # Add regression equation, R-squared, and p-value as text annotation
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")),
    size = 2.5
  )


# --- Figure S3b: Residence Time vs. Isolation (Invaded sites only) ---

# Load network data and filter to include only invaded sites (invasion == "Y")
plot_dat <- read_csv("network_data.csv") %>%
  filter(invasion == "Y")

# Apply log10 transformation to the isolation variable
plot_dat$isolation <- log10(plot_dat$isolation)

# Create a scatter plot with a linear regression line for invaded sites
figs2b <- ggplot(data = plot_dat, aes(x = isolation, y = residence_time)) +
  # Scatter points: all black, semi-transparent
  geom_point(size = 3, shape = 21, alpha = 5/10, fill = "black", colour = "black") +
  # Axis labels
  labs(y = "Residence time since introduction (year)", x = "Distance to closest mainland [log10(km)]") +
  # Add a linear regression line (method = "lm")
  stat_smooth(method = "lm", size = 1.5, linetype = 1, colour = "black") +
  # Add regression equation, R-squared, and p-value as text annotation
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~~")),
    size = 2.5
  )


# --- Assemble and Save the Combined Figure ---

# Combine figs2a and figs2b into a single figure with 2 columns
# Use "AUTO" to automatically generate labels "A" and "B"
figs2 <- plot_grid(figs2a, figs2b, ncol = 2, labels = "AUTO")

# Display the combined figure
figs2

# Save the combined figure to a file
ggsave(
  filename = paste0('figure/figure_s2.png'),
  plot = figs2,
  width = 6,  # Width in inches
  height = 3  # Height in inches
)
