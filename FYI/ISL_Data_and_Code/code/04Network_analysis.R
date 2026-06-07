# Clear the entire workspace to avoid conflicts with previous variables
rm(list = ls())

# Load required R packages
library(bipartite)  # For bipartite network analysis
library(tidyverse)  # For data manipulation and visualization (includes dplyr, ggplot2, etc.)
library(reshape2)   # For reshaping data from long to wide format

# Load the site-specific network data from a CSV file
# The data is in long format with columns: prey, species, RRA (relative read abundance)
site_network <- read_csv("data/example_site1_JT_site060.csv")

# Reshape the data from long to wide format to create the interaction matrix
# - 'prey' as rows, 'species' (predators) as columns
# - 'RRA' (relative read abundance) as the values
# - 'fun.aggregate = sum' sums RRA values for any duplicate prey-predator pairs
network_matrix <- site_network %>%
  reshape2::acast(formula = prey ~ species, fun.aggregate = sum, value.var = "RRA") %>%
  as.matrix()

# --- Null Model Analysis ---

# The bipartite::nullmodel function requires integer counts.
# Therefore, multiply the matrix by 100 and round to convert proportions to integers.
site_matrix2 <- round(network_matrix * 100)

# Generate 1000 random null model networks using 'method = 1' (r2dtable)
# This method preserves the row and column totals (marginal sums) of the original matrix,
# randomly redistributing the interactions while maintaining species' total interaction frequencies.
nulls <- nullmodel(site_matrix2, method = 1)

# Convert the null model matrices back to proportions by dividing by 100,
# to match the scale of the original 'network_matrix'.
for (i in 1:1000) {
  nulls[[i]] <- nulls[[i]] / 100
}

# --- Network-Level Analysis ---

# Calculate the observed connectance of the network.
# Connectance is the proportion of realized interactions out of all possible interactions.
net_connect <- networklevel(network_matrix, index = "connectance")

# Calculate connectance for each of the 1000 null model networks.
Inulls_connect <- sapply(nulls, function(x) networklevel(x, index = "connectance"))

# Compute the z-score for the observed connectance.
# The z-score indicates how many standard deviations the observed value is from the mean of the null distribution.
z_connect <- (net_connect - mean(Inulls_connect)) / sd(Inulls_connect)

# Calculate the two-tailed p-value for the z-score.
# This p-value tests whether the observed connectance is significantly different from the null expectation.
p_connect <- 2 * pnorm(-abs(z_connect))

# Calculate the observed modularity (Q) of the network.
# Modularity measures the extent to which the network is divided into densely connected modules.
# 'computeModules' uses an algorithm to find the optimal modular structure.
Modularity_Q <- c(Modularity.Q = computeModules(network_matrix)@likelihood)

# Calculate modularity for each of the 1000 null model networks.
Inulls_Q <- sapply(nulls, function(x) computeModules(x)@likelihood)

# Compute the z-score for the observed modularity.
z_Modularity_Q <- (Modularity_Q - mean(Inulls_Q)) / sd(Inulls_Q)

# Calculate the two-tailed p-value for the modularity z-score.
p_Q <- 2 * pnorm(-abs(z_Modularity_Q))

# Calculate the observed weighted nestedness (WNODF) of the network.
# Nestedness measures the extent to which specialist species interact with subsets of the partners of generalist species.
# 'weighted = T' accounts for interaction strengths (RRA values).
WNODF <- nest.smdm(network_matrix, weighted = T)$WNODFmatrix

# Calculate weighted nestedness for each of the 1000 null model networks.
Inulls_WNODF <- sapply(nulls, function(x) nest.smdm(x, weighted = T)$WNODFmatrix)

# Compute the z-score for the observed nestedness.
z_WNODF <- (WNODF - mean(Inulls_WNODF)) / sd(Inulls_WNODF)

# Calculate the two-tailed p-value for the nestedness z-score.
p_WNODF <- 2 * pnorm(-abs(z_WNODF))


# --- Species-Level (Centrality) Analysis ---

# Calculate the observed 'degree' for each predator species (level = "higher").
# Degree is the number of different prey species each predator interacts with.
Predator_degree <- specieslevel(network_matrix, index = c("degree"), level = "higher") %>%
  as.data.frame()

# Calculate degree for each predator in each of the null model networks.
Inulls_degree <- sapply(nulls, function(x) specieslevel(x, index = c("degree"), level = "higher"))

# For each predator, compute the mean degree across all null model networks.
Predator_degree$mean_null_degree <- rowMeans(data.frame(Inulls_degree))

# For each predator, compute the standard deviation of degree across all null model networks.
Predator_degree$sd_null_degree <- apply(data.frame(Inulls_degree), 1, sd)

# Compute the z-score for each predator's observed degree.
# This z-score indicates whether a predator has significantly more or fewer interactions than expected by chance.
Predator_degree$z_degree <- (Predator_degree$degree - Predator_degree$mean_null_degree) / Predator_degree$sd_null_degree


# Calculate the observed 'weighted closeness' for each predator species.
# Closeness centrality measures how quickly a node can reach all other nodes in the network.
# The weighted version accounts for interaction strengths.
Predator_closeness <- specieslevel(as.matrix(network_matrix), index = c("closeness"), level = "higher") %>%
  as.data.frame() %>%
  dplyr::select(weighted.closeness) # Extract only the weighted closeness column

# Calculate weighted closeness for each predator in each of the null model networks.
Inulls_closeness <- sapply(nulls, function(x) specieslevel(x, index = c("closeness"), level = "higher"))

# The result is a list of data frames. Extract the 'weighted.closeness' values for each null model.
Inulls_closeness_list <- lapply(1:1000, function(i) Inulls_closeness[, i]$weighted.closeness)

# For each predator, compute the mean weighted closeness across all null model networks.
Predator_closeness$mean_null_closeness <- rowMeans(data.frame(Inulls_closeness_list))

# For each predator, compute the standard deviation of weighted closeness across all null model networks.
Predator_closeness$sd_null_closeness <- apply(data.frame(Inulls_closeness_list), 1, sd)

# Compute the z-score for each predator's observed weighted closeness.
# This indicates whether a predator is significantly more or less centrally located (in terms of closeness) than expected by chance.
Predator_closeness$z_closeness <- (Predator_closeness$weighted.closeness - Predator_closeness$mean_null_closeness) / Predator_closeness$sd_null_closeness
```