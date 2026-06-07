# Load required R packages
library(piecewiseSEM)  # For fitting piecewise structural equation models (SEM)
library(tidyverse)     # For data manipulation (dplyr) and visualization (ggplot2)
library(inspectdf)     # For inspecting data frames (not explicitly used in this script, but loaded)
library(DHARMa)        # For model diagnostics (not explicitly used here, but loaded)
library(glmmTMB)       # For fitting generalized linear mixed-effects models (GLMMs), used within SEM

# --- Load and Prepare Data ---

# Load the network data containing metrics and predictor variables
network_dat <- read_csv("data/network_data.csv")

# Define a custom function "scale2" for z-transformation (standardization) of continuous variables
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# Create a new dataframe "predictors_std" to store standardized predictors
predictors_std <- network_dat

# Standardize each continuous predictor (apply log transformations first if needed to normalize distribution)
predictors_std$area <- scale2(log10(predictors_std$area))  # Area: log10 transformation + z-score
predictors_std$isolation <- scale2(log10(predictors_std$isolation))  # Isolation: log10 + z-score
predictors_std$residence_time <- scale2(predictors_std$residence_time)  # Residence time: direct z-score
predictors_std$bullfrog_abundance <- scale2(log(predictors_std$bullfrog_abundance + 1))  # Bullfrog abundance: log(x+1) (avoid 0) + z-score
predictors_std$native_abundance <- scale2(log(predictors_std$native_abundance + 1))  # Native species abundance: log(x+1) + z-score
predictors_std$crayfish_abundance <- scale2(log(predictors_std$crayfish_abundance + 1))  # Crayfish abundance: log(x+1) + z-score
predictors_std$water_area <- scale2(log10(predictors_std$water_area))  # Water area: log10 + z-score
predictors_std$human_activity <- scale2(predictors_std$human_activity)  # Human activity: direct z-score
predictors_std$vegetation_coverage <- scale2(predictors_std$vegetation_coverage)  # Vegetation coverage: direct z-score

# Convert categorical variables to factor type (required for mixed-effects model random/fixed effects)
predictors_std$island <- as.factor(predictors_std$island)  # Island ID (categorical)
predictors_std$year <- as.factor(predictors_std$year)  # Year (categorical)

# --- Add Spatial Correlation Component ---

# Create a spatial factor variable "pos" from x and y coordinates
# This is used to model spatial autocorrelation within islands using glmmTMB's covariance structures
predictors_std$pos = numFactor(predictors_std$x, predictors_std$y) 
# Reference: https://glmmtmb.github.io/glmmTMB/articles/covstruct.html#spatial-correlations

# --- Define Optimizer Controls for glmmTMB ---

# Define different optimizers to improve model convergence, especially for complex models
opti_BFGS <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
opti_CG <- glmmTMBControl(optimizer = optim, optArgs = list(method = "CG"))
opti_nlinb <- glmmTMBControl(optimizer = nlminb, optCtrl = list(iter.max = 1e4, eval.max = 1e4))


# --- Structural Equation Model (SEM) 1: Connectance ---

# Set NA handling to omit missing values
options(na.action = "na.omit")

# Define the first piecewise SEM using 'psem'
# The model consists of two regression equations (GLMMs) representing the hypothesized causal relationships:
# 1. Connectance is affected by bullfrog abundance and other predictors, with random and spatial effects.
# 2. Bullfrog abundance is affected by residence time and other predictors, with random effects.
mod1 <- psem(
  # Submodel 1: Response = Connectance (endogenous variable)
  glmmTMB(Connectance ~ bullfrog_abundance + residence_time + native_abundance + area + isolation
          + water_area + vegetation_coverage + human_activity + crayfish_abundance
          + (1|island) + (1|year) 
          + exp(pos + 0|island),  # Exponential spatial autocorrelation structure within islands
          family = gaussian(), data = predictors_std),
  
  # Submodel 2: Response = bullfrog_abundance (endogenous variable predicting Connectance)
  glmmTMB(bullfrog_abundance ~ residence_time + native_abundance + area + isolation + water_area
          + (1|island) + (1|year),
          data = predictors_std, control = opti_BFGS)
)

# --- Model 1 Summary and Diagnostics ---

# Print a detailed summary of the SEM, including standardized path coefficients, p-values, and model fit statistics
summary(mod1, .progressBar = T) 

# Perform a Fisher's C test to assess the overall goodness-of-fit of the SEM.
# A non-significant p-value (p > 0.05) indicates that the model fits the data well (i.e., the data are consistent with the model structure).
fisherC(mod1) 

# Perform d-separation tests to evaluate the conditional independence claims of the model.
# Non-significant tests support the model's assumptions.
dSep(mod1)

# Generate a basic visual plot of the SEM path diagram (for quick inspection)
plot(mod1)

# --- Extract and Format Results for Model 1 ---

# Extract coefficients from the SEM and format the dataframe for export
sem_dt1 <- coefs(mod1) %>%
  select(-9) %>%  # Remove the 9th column (typically the 'Std.Est.2.5%' or similar)
  filter(
    # Filter out rows representing residual covariances (denoted by "~~")
    str_detect(Response, "~~", negate = TRUE) 
  ) %>%
  # Reorder columns for clarity and ease of use in visualization software
  dplyr::relocate(
    from = Predictor,
    to = Response,
    weight = Std.Estimate,
    p = P.Value
  )


# --- Structural Equation Model (SEM) 2: Modularity ---

# Define the second piecewise SEM, identical in structure to mod1 but with Modularity as the response variable
mod2 <- psem(
  glmmTMB(Modularity ~ bullfrog_abundance + residence_time + native_abundance + area + isolation
          + water_area + vegetation_coverage + human_activity + crayfish_abundance
          + (1|island) + (1|year) 
          + exp(pos + 0|island),
          family = gaussian(), data = predictors_std),
  
  glmmTMB(bullfrog_abundance ~ residence_time + native_abundance + area + isolation + water_area
          + (1|island) + (1|year),
          data = predictors_std, control = opti_BFGS)
)

# --- Model 2 Summary and Diagnostics ---
summary(mod2, .progressBar = T) 
fisherC(mod2) 
dSep(mod2)
plot(mod2)

# --- Extract and Format Results for Model 2 ---
sem_dt2 <- coefs(mod2) |>
  select(-9) %>%
  filter(
    str_detect(Response, "~~", negate = TRUE) 
  ) %>%
  dplyr::relocate(
    from = Predictor,
    to = Response,
    weight = Std.Estimate,
    p = P.Value
  )


# --- Structural Equation Model (SEM) 3: Nestedness ---

# Define the third piecewise SEM, with Nestedness as the response variable
# Note: `ziformula=~0` is redundant for a Gaussian family and can be omitted.
mod3 <- psem(
  glmmTMB(Nestedness ~ bullfrog_abundance + residence_time + native_abundance + area + isolation
          + water_area + vegetation_coverage + human_activity + crayfish_abundance
          + (1|island) + (1|year) 
          + exp(pos + 0|island),
          ziformula = ~0, family = gaussian(), data = predictors_std),
  
  glmmTMB(bullfrog_abundance ~ residence_time + native_abundance + area + isolation + water_area
          + (1|island) + (1|year),
          data = predictors_std, control = opti_BFGS)
)

# --- Model 3 Summary and Diagnostics ---
summary(mod3, .progressBar = T) 
fisherC(mod3) 
dSep(mod3)
plot(mod3)

# --- Extract and Format Results for Model 3 ---
sem_dt3 <- coefs(mod3) %>%
  select(-9) %>%
  filter(
    str_detect(Response, "~~", negate = TRUE) 
  ) %>%
  dplyr::relocate(
    from = Predictor,
    to = Response,
    weight = Std.Estimate,
    p = P.Value
  )


# --- Combine and Export Results for Visualization ---

# The following code prepares a comprehensive table for creating a publication-ready
# SEM path diagram in external software like Adobe Illustrator.

# Define custom orders for grouping and arranging variables in the visualization
custom_order1 <- c("Invasion characteristics", "Island characteristics", "Native characteristics", "Habitat characteristics")
custom_order2 <- c("bullfrog_abundance", "Connectance", "Modularity", "Nestedness")
custom_order3 <- c("residence_time", "bullfrog_abundance", "isolation", "area",
                   "native_abundance", "water_area", "vegetation_coverage",
                   "human_activity", "crayfish_abundance")

# Combine results from all three SEMs into a single dataframe
dt <- rbind(sem_dt1, sem_dt2, sem_dt3) %>%
  # Create a 'Groups' column to categorize predictors for better visualization
  mutate(Groups = case_when(
    from == "residence_time" ~ "Invasion characteristics",
    from == "bullfrog_abundance" ~ "Invasion characteristics",
    from == "isolation" ~ "Island characteristics",
    from == "area" ~ "Island characteristics",
    from == "native_abundance" ~ "Native characteristics",
    from == "water_area" ~ "Habitat characteristics",
    from == "vegetation_coverage" ~ "Habitat characteristics",
    from == "human_activity" ~ "Habitat characteristics",
    from == "crayfish_abundance" ~ "Habitat characteristics"
  )) %>%
  # Select and reorder columns for clarity
  select(Groups, from, to, weight, p, Estimate, DF, Std.Error) %>%
  distinct() # Remove any duplicate rows

# Apply the custom factor orders to ensure variables appear in the desired order in plots/software
dt$Groups <- factor(dt$Groups, levels = custom_order1)
dt$to <- factor(dt$to, levels = custom_order2)
dt$from <- factor(dt$from, levels = custom_order3)

# Sort the dataframe by the response and then predictor variables
dt <- dt[order(dt$to), ]
dt <- dt[order(dt$from), ]

# Format variable names for readability (convert to title case and replace underscores with spaces)
dt$from <- str_to_title(dt$from)
dt$from <- str_replace(dt$from, "_", " ")

# Rename columns to more descriptive names for the final output table
colnames(dt) <- c("Groups", "Predictor", "Response", "Standard estimate", "P value", "Parameter estimate", "DF", "Standard error")

# Export the formatted table to a CSV file for use in external visualization software
write_excel_csv(dt, "data/sem_sum_table.csv")
