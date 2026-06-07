```r
# Load required R packages for data processing, visualization, modeling, and diagnostics
library(tidyverse)    # For data cleaning, manipulation, and visualization (core package collection)
library(ggpubr)       # For creating publication-ready ggplot2 plots
library(ggpmisc)      # For adding statistical annotations (e.g., equations, R²) to ggplot2 plots
library(cowplot)      # For combining multiple ggplot2 plots into a single figure
library(car)          # For variance analysis (ANOVA) and regression diagnostics
library(glmmTMB)      # For fitting generalized linear mixed-effects models (GLMMs) with TMB backend
library(DHARMa)       # For residual diagnostics of hierarchical (mixed) models

# Read the network data from CSV file (path: "data/network_data.csv")
dat <- read_csv("data/network_data.csv")

################################################################################
### Section 1: Check for Spatial Autocorrelation
################################################################################
library(sp)           # For handling spatial data structures (e.g., spatial points)
library(spdep)        # For spatial dependence analysis (e.g., Moran's I test)
library(ggplot2)      # For basic data visualization (redundant with tidyverse but kept for clarity)
library(sf)           # For modern spatial data formats and operations (simple features)

# Convert the dataframe to a spatial point object (sp package format) using x (longitude) and y (latitude)
coordinates(dat) <- ~x + y  
# Define the coordinate reference system (CRS) as WGS84 (global geographic coordinate system)
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")  

# Convert geographic coordinates (WGS84) to UTM coordinates (meters) to ensure correct distance units
df_utm <- spTransform(dat, CRS("+proj=utm +zone=51 +datum=WGS84"))  

# Build a neighbor list using K-nearest neighbors (K=5: each point is connected to its 5 closest points)
nb <- knn2nb(knearneigh(df_utm, k = 5))  
# Convert the neighbor list to a spatial weight matrix (style = "W": row-standardized weights)
weight_matrix <- nb2listw(nb, style = "W")  

# Perform Moran's I test for spatial autocorrelation (response variable: Connectance)
moran_result <- moran.test(dat$Connectance, weight_matrix, na.action = na.omit)
print(moran_result)  # Result note: p-value = 0.066

# Perform Moran's I test for spatial autocorrelation (response variable: Modularity)
moran_result <- moran.test(dat$Modularity, weight_matrix)
print(moran_result)  # Result note: p-value = 0.030

# Perform Moran's I test for spatial autocorrelation (response variable: Nestedness)
moran_result <- moran.test(dat$Nestedness, weight_matrix)
print(moran_result)  # Result note: p-value = 0.006

# Conclusion on spatial autocorrelation:
# Spatial autocorrelation exists: p < 0.05 and Moran's I > 0 indicate positive autocorrelation


################################################################################
### Section 2: Data Normalization Transformations
################################################################################
# Re-read the original data (to avoid spatial object formatting from Section 1)
dat <- read_csv("data/network_data.csv")

# Define a custom function "scale2" for z-transformation (standardization) of continuous variables
# z-transformation: (value - mean) / standard deviation; handles NA values if specified
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# Create a new dataframe "predictors_std" for standardized predictors
predictors_std <- dat

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

################################################################################
### Section 3: Add Spatial Correlation Variable
################################################################################
# Create a spatial factor variable "pos" using x (longitude) and y (latitude) coordinates
# Used to model spatial autocorrelation in glmmTMB (see reference for spatial covariance structure)
predictors_std$pos = numFactor(predictors_std$x, predictors_std$y) 
# Reference: glmmTMB documentation on spatial covariance structures
# https://glmmtmb.github.io/glmmTMB/articles/covstruct.html#spatial-correlations


################################################################################
### Section 4: Define Optimizer Controls for glmmTMB Models
################################################################################
# Define model control parameters with different optimizers (to improve model convergence)
# 1. BFGS optimizer (commonly used for smooth objective functions)
opti_BFGS <- glmmTMBControl(optimizer = optim,
                            optArgs = list(method = "BFGS"))

# 2. Conjugate Gradient (CG) optimizer (effective for high-dimensional data)
opti_CG <- glmmTMBControl(optimizer = optim,
                          optArgs = list(method = "CG"))

# 3. nlminb optimizer (non-linear minimization with bounds; increased iteration limits)
opti_nlinb <- glmmTMBControl(optimizer = nlminb,
                             optCtrl = list(iter.max = 1e4, eval.max = 1e4))  # 10,000 iterations/evaluations

# Response variables for the three main models (network structure metrics)
# - Connectance: Proportion of potential links realized in the network
# - Modularity: Degree of network division into densely connected submodules
# - Nestedness: Degree of specialist species interactions being subsets of generalist interactions


################################################################################
### Section 5: Check Multicollinearity Among Predictors
################################################################################
# Multicollinearity check for Connectance (response variable)
# Fit linear model with all continuous predictors
lmt <- lm(Connectance ~ bullfrog_abundance + residence_time + native_abundance
          + area + isolation + water_area + vegetation_coverage + human_activity 
          + crayfish_abundance, data = predictors_std %>% filter(!is.na(Connectance)))
# Check Variance Inflation Factor (VIF): VIF < 4 indicates no severe multicollinearity
performance::check_collinearity(lmt)  # Result: All VIF < 4

# Multicollinearity check for Modularity (response variable)
lmt <- lm(Modularity ~ bullfrog_abundance + residence_time + native_abundance
          + area + isolation + water_area + vegetation_coverage + human_activity 
          + crayfish_abundance, data = predictors_std)
performance::check_collinearity(lmt)  # Result: All VIF < 4

# Multicollinearity check for Nestedness (response variable)
lmt <- lm(Nestedness ~ bullfrog_abundance + residence_time + native_abundance
          + area + isolation + water_area + vegetation_coverage + human_activity 
          + crayfish_abundance, data = predictors_std)
performance::check_collinearity(lmt)  # Result: All VIF < 4


################################################################################
### Section 6: Model 1 - Connectance (Response Variable)
################################################################################
# Set NA handling rule: Fail if NA values are present (ensures data quality)
options(na.action = "na.fail")

# Subset data to exclude NA values in Connectance (response variable)
predictors_std1 <- predictors_std %>% filter(!is.na(Connectance))

# Fit mixed-effects model for Connectance
# - Fixed effects: residence_time, native_abundance, area×bullfrog_abundance (interaction), 
#   isolation×bullfrog_abundance (interaction), water_area, vegetation_coverage, human_activity, crayfish_abundance
# - Random effects: (1|island) (island-specific intercept), (1|year) (year-specific intercept)
# - Spatial structure: exp(pos + 0 | island) (exponential spatial covariance within islands)
# - Distribution: Gaussian (normal distribution for continuous response)
# - Optimizer: BFGS (defined in Section 4)
Connectance_model <- glmmTMB(Connectance ~ residence_time 
                             + native_abundance
                             + (area + isolation)*bullfrog_abundance  # Interaction terms: area×bullfrog, isolation×bullfrog
                             + water_area 
                             + vegetation_coverage
                             + human_activity 
                             + crayfish_abundance
                             + (1|island)  # Random intercept for island
                             + (1|year)    # Random intercept for year
                             + exp(pos + 0 | island)  # Exponential spatial autocorrelation within islands
                             ,
                             family = gaussian(), 
                             data = predictors_std1,
                             control = opti_BFGS)

# Model summary
summary(Connectance_model)  # Print detailed model summary (coefficients, SE, p-values, random effects)

################################################################################
### Section 7: Model 1 Diagnostics (Connectance Model)
################################################################################
# Simulate residuals for model diagnostics (DHARMa package)
sim_res <- simulateResiduals(Connectance_model)
plot(sim_res)  # Plot residual diagnostics (QQ plot, residual vs. predicted, etc.)

# Test for overdispersion (common in count data; less critical for Gaussian models)
testDispersion(Connectance_model)  

# Test for zero-inflation (whether observed zeros exceed model predictions)
# Note: Significant result (p < 0.05) indicates zero-inflation: model predicts fewer zeros than observed
testZeroInflation(Connectance_model)  

# Test for spatial autocorrelation in residuals (verify if spatial structure is captured)
testSpatialAutocorrelation(sim_res, x = predictors_std1$x, y = predictors_std1$y)  # Result: p-value = 0.1238 (no residual spatial autocorrelation)


################################################################################
### Section 8: Model Selection and Averaging (Connectance)
################################################################################
library(MuMIn)  # For model selection and averaging

# Perform model dredging (all possible subsets of fixed effects)
Connectance_model_dredge <- MuMIn::dredge(Connectance_model)
# Remove models with NA delta values (poorly fitted models)
Connectance_model_dredge <- Connectance_model_dredge[!is.na(Connectance_model_dredge$delta), ]

# Select top models with cumulative AIC weight ≤ 0.95 (captures 95% of model uncertainty)
Connectance_model_dredge_w95 <- get.models(Connectance_model_dredge, subset = cumsum(weight) <= .95)
# Perform model averaging on selected top models
Connectance_model_dredge.avgmod.top5 <- model.avg(Connectance_model_dredge_w95)

# Model averaging output
summary(Connectance_model_dredge.avgmod_95)  # Print averaged model summary
confint(Connectance_model_dredge.avgmod_95)  # 95% confidence intervals for averaged coefficients

# Extract model weights (sw = variable importance weight) and format
con_avg_sw <- data.frame(sw(Connectance_model_dredge.avgmod_95)) %>%  # Note: Typo in original ("avgmod_95" vs "avgmod.top5")
  mutate(
    variables = str_to_title(c(str_replace(str_replace(rownames(.), "cond[(]", ""), "[)]", ""))),  # Clean variable names
    network = "Connectance"  # Label network metric
  )
colnames(con_avg_sw) <- c("sw", "variables", "network")  # Rename columns

# Extract averaged coefficients, 95% CI, p-values, and merge with weights
dd_con_avg_cof <- cbind(
  confint(Connectance_model_dredge.avgmod_95),  # 95% CI (lower, upper)
  coef(Connectance_model_dredge.avgmod_95),     # Averaged standardized coefficients
  summary(Connectance_model_dredge.avgmod_95)$coefmat.subset[, 5]  # p-values
)[-1, ] %>%  # Remove intercept row
  as.data.frame() %>%
  rename(stand_cof = 3, p = 4) %>%  # Rename coefficient and p-value columns
  mutate(
    variables = str_to_title(c(str_replace(str_replace(rownames(.), "cond[(]", ""), "[)]", ""))),  # Clean variable names
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rownames_to_column(var = "row.names") %>%  # Add row names as a column
  left_join(con_avg_sw)  # Merge with variable weights

# Set factor levels for variables (for consistent ordering in plots)
dd_con_avg_cof$variables <- factor(dd_con_avg_cof$variables, 
                                  levels = rev(c("Isolation", "Area", "Residence_time", "Bullfrog_abundance",
                                                 "Native_abundance", "Water_area", "Vegetation_coverage",
                                                 "Human_activity", "Crayfish_abundance",
                                                 "Bullfrog_abundance:Isolation", "Area:Bullfrog_abundance")))


################################################################################
### Section 9: Model 2 - Modularity (Response Variable)
################################################################################
# Fit mixed-effects model for Modularity
# - Fixed effects: Same as Connectance model (excludes bullfrog_abundance main effect, per original code)
# - Random effects and spatial structure: Same as Connectance model
# - Distribution: Gaussian
Modularity_model <- glmmTMB(Modularity ~ 
                              residence_time
                            + native_abundance
                            + (area + isolation)*bullfrog_abundance  # Interactions: area×bullfrog, isolation×bullfrog
                            + water_area 
                            + vegetation_coverage
                            + human_activity 
                            + crayfish_abundance
                            + (1|island)  # Random intercept for island
                            + (1|year)    # Random intercept for year
                            + exp(pos + 0 | island),  # Exponential spatial autocorrelation
                            family = gaussian(),
                            data = predictors_std)  # No NA filtering (original code)

# Model summary
summary(Modularity_model)


################################################################################
### Section 10: Model 2 Diagnostics (Modularity Model)
################################################################################
sim_res <- simulateResiduals(Modularity_model)
plot(sim_res)  # Overall residual diagnostics
testDispersion(Modularity_model)  # Test overdispersion
# Test zero-inflation: p < 0.05 indicates zero-inflation (model predicts fewer zeros than observed)
testZeroInflation(Modularity_model)  
# Test residual spatial autocorrelation
testSpatialAutocorrelation(sim_res, x = predictors_std$x, y = predictors_std$y)


################################################################################
### Section 11: Model Selection and Averaging (Modularity)
################################################################################
options(na.action = "na.fail")
library(MuMIn)

# Model dredging and filtering
Modularity_model_dredge <- MuMIn::dredge(Modularity_model)
Modularity_model_dredge$delta  # View delta AIC values
Modularity_model_dredge <- Modularity_model_dredge[!is.na(Modularity_model_dredge$delta), ]

# Select top models (cumulative weight ≤ 0.95) and average
Modularity_model_dredge_w95 <- get.models(Modularity_model_dredge, subset = cumsum(weight) <= .95)
Modularity_model_dredge.avgmod_95 <- model.avg(Modularity_model_dredge_w95)

# Averaged model output
summary(Modularity_model_dredge.avgmod_95)
confint(Modularity_model_dredge.avgmod_95)

# Extract variable weights and format
mod_avg_sw <- data.frame(sw(Modularity_model_dredge.avgmod_95)) %>%
  mutate(
    variables = str_to_title(c(str_replace(str_replace(rownames(.), "cond[(]", ""), "[)]", ""))),
    network = "Modularity"
  )
colnames(mod_avg_sw) <- c("sw", "variables", "network")

# Extract averaged coefficients, CI, p-values, and merge with weights
dd_mod_avg_cof <- cbind(
  confint(Modularity_model_dredge.avgmod_95),
  coef(Modularity_model_dredge.avgmod_95),
  summary(Modularity_model_dredge.avgmod_95)$coefmat.subset[, 5]
)[-1, ] %>%
  as.data.frame() %>%
  rename(stand_cof = 3, p = 4) %>%
  mutate(
    variables = str_to_title(c(str_replace(str_replace(rownames(.), "cond[(]", ""), "[)]", ""))),
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rownames_to_column(var = "row.names") %>%
  left_join(mod_avg_sw)

# Set variable factor levels for plotting
dd_mod_avg_cof$variables <- factor(dd_mod_avg_cof$variables, 
                                  levels = rev(c("Isolation", "Area", "Residence_time", "Bullfrog_abundance",
                                                 "Native_abundance", "Water_area", "Vegetation_coverage",
                                                 "Human_activity", "Crayfish_abundance",
                                                 "Bullfrog_abundance:Isolation", "Area:Bullfrog_abundance")))



################################################################################
### Section 12: Model 3 - Nestedness (Response Variable)
################################################################################
# Fit mixed-effects model for Nestedness
# - Fixed effects: Same as previous models (explicitly lists interaction terms)
# - Random effects and spatial structure: Same as previous models
# - Distribution: Gaussian; Optimizer: nlminb (for better convergence)
# - REML = TRUE: Restricted Maximum Likelihood (reduces bias in variance estimates)
Nestedness_model <- glmmTMB(Nestedness ~ 
                              residence_time
                            + native_abundance
                            + area + isolation + bullfrog_abundance  # Main effects
                            + area:bullfrog_abundance  # Explicit interaction: area×bullfrog
                            + isolation:bullfrog_abundance  # Explicit interaction: isolation×bullfrog
                            + water_area 
                            + vegetation_coverage
                            + human_activity 
                            + crayfish_abundance
                            + (1|island)  # Random intercept for island
                            + (1|year)    # Random intercept for year
                            + exp(pos + 0 | island),  # Exponential spatial autocorrelation
                            REML = TRUE,
                            family = gaussian(),
                            data = predictors_std,
                            control = opti_nlinb)  # Use nlminb optimizer

# Model summary
summary(Nestedness_model)


################################################################################
### Section 13: Model 3 Diagnostics (Nestedness Model)
################################################################################
sim_res <- simulateResiduals(Nestedness_model)
plot(sim_res)  # Overall residual diagnostics
testDispersion(Nestedness_model)  # Test overdispersion
# Test zero-inflation: p < 0.05 indicates zero-inflation
testZeroInflation(Nestedness_model)  
# Test residual spatial autocorrelation
testSpatialAutocorrelation(sim_res, x = predictors_std$x, y = predictors_std$y)


################################################################################
### Section 14: Model Selection and Averaging (Nestedness)
################################################################################
library(MuMIn)

# Model dredging and filtering
Nestedness_model_dredge <- MuMIn::dredge(Nestedness_model)
Nestedness_model_dredge <- Nestedness_model_dredge[!is.na(Nestedness_model_dredge$delta), ]

# Select top models (cumulative weight ≤ 0.95) and average
Nestedness_model_dredge_w95 <- get.models(Nestedness_model_dredge, subset = cumsum(weight) <= .95)
Nestedness_model_dredge.avgmod_95 <- model.avg(Nestedness_model_dredge_w95)

# Averaged model output
summary(Nestedness_model_dredge.avgmod_95)
confint(Nestedness_model_dredge.avgmod_95)

# Extract variable weights and format
nest_avg_sw <- data.frame(sw(Nestedness_model_dredge.avgmod_95)) %>%
  mutate(
    variables = str_to_title(c(str_replace(str_replace(rownames(.), "cond[(]", ""), "[)]", ""))),
    network = "Nestedness"
  )
colnames(nest_avg_sw) <- c("sw", "variables", "network")

# Extract averaged coefficients, CI, p-values, and merge with weights
dd_nest_avg_cof <- cbind(
  confint(Nestedness_model_dredge.avgmod_95),
  coef(Nestedness_model_dredge.avgmod_95),
  summary(Nestedness_model_dredge.avgmod_95)$coefmat.subset[, 5]
)[-1, ] %>%
  as.data.frame() %>%
  rename(stand_cof = 3, p = 4) %>%
  mutate(
    variables = str_to_title(c(str_replace(str_replace(rownames(.), "cond[(]", ""), "[)]", ""))),
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rownames_to_column(var = "row.names") %>%
  left_join(nest_avg_sw)

# Set variable factor levels for plotting
dd_nest_avg_cof$variables <- factor(dd_nest_avg_cof$variables, 
                                   levels = rev(c("Isolation", "Area", "Residence_time", "Bullfrog_abundance",
                                                  "Native_abundance", "Water_area", "Vegetation_coverage",
                                                  "Human_activity", "Crayfish_abundance",
                                                  "Bullfrog_abundance:Isolation", "Area:Bullfrog_abundance")))



################################################################################
### Section 15: Combine Results and Export
################################################################################
# Merge averaged coefficients from all three models (note: original typo "dd_nest_con_cof" → corrected to "dd_con_avg_cof")
dt <- rbind(dd_con_avg_cof, dd_mod_avg_cof, dd_nest_avg_cof)

# Export combined results to CSV file (for Figure 5: effect size plots)
write_excel_csv(dt, "data/fig4-glmm_effect_size.csv")

```
