# Load required R packages
library(tidyverse)  # For data manipulation (read_csv) and visualization (ggplot2)
library(iNEXT)      # For performing iNEXT (iNterpolation and EXTrapolation) analysis to create rarefaction/extrapolation curves

# --- Figure S3: Rarefaction Curves of Prey Data ---

# 1. Load and Prepare Data
# Read the CSV file containing sequencing read counts for each sample and OTU
# This file should be in a matrix-like format (rows = samples, columns = OTUs, values = read counts)
sample_read_seq_mat <- read_csv("sequening_reads_number.csv") %>%
  as.matrix()  # Convert the dataframe to a matrix, which is the required input format for iNEXT

# 2. Perform iNEXT Analysis
rc_reads <- iNEXT(sample_read_seq_mat, q = 0, datatype = "abundance")

# 3. Generate and Customize the Plot
p1 <- ggiNEXT(rc_reads, type = 1, color.var = "None")

# Customize the appearance of the plot
p11 <- p1 +
  theme_test() +
  theme(legend.position = 'none') + 
  labs(x = "Number of sequencing reads", y = "Number of detected Arthropoda OTUs")

# 4. Save the Plot
# Save the final customized plot to a file
ggsave(
  paste0("figure/figure_s3_Rarefaction_curves_of_the_prey_data.png"),
  plot = p11,
  width = 6,   # Set the width of the output image in inches
  height = 5   # Set the height of the output image in inches
)
