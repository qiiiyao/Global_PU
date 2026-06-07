This repository hosts the main code for generating the results for "Drought legacies delay spring green up in northern ecosystems"

1. System requirements
operating systems: windows 10
software version: Matlab 2022b, Python 3.10.9., R 4.1.3

2. Demo
instructions to run on data
    1. extracting plant phenology for the GIMMS NDVI
    2. identifying drought events for the GIMMS NDVI and SM
    3. analyzing drought legacy effect on spring phenology 
    4. exploring the underlying mechanisms of drought legacy using Random Forest (RF) and Structural equation modeling (SEM)

Expected output
     1. plant phenology 
     2. drought types
     3. drought-induced SOS changes
     4. RF variable importance and partial dependent plot; standardized path coefficients of path diagrams

3. Instructions for use
    1. imput  GIMMS NDVI time series
    2. imput  monthly SM, NDVI and yearly SOS and EOS
    3. imput  observated SOS and modeled SOS
    4. imput  RF variable and SEM variable using RF_prepare.m and SEM_prepare.m