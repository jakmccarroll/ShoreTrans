FILES IN THE ZENODO REPOSITORY: https://doi.org/10.5281/zenodo.15161568
- This repository contains codes and data for the article under review: "Morphologically adaptive modelling of sea level rise induced coastal erosion impact for south-east Australia"


QUICK INSTRUCTIONS
1. Download and unzip "_tutorial_scripts+data_22Jun2025.zip"
2. Download and unzip "_ShoreTrans_functions_22Jun2025.zip" ---> add all functions to the Matlab path
3. Download and unzip "_profile_prep+classify_functions_22Jun2025.zip" ---> add all functions to the Matlab path
3. Open script "../tutorial/scripts/tut01_prep_classify.m" --> read through instructions, run the profile classification for the 7 example profiles
4. Open script "../tutorial/scripts/tut02_ST_run.m" --> run ShoreTrans for the 7 example profiles
5. To examine model outputs (cross-sections) for any profile from the full Victorian ShoreTrans model, download all data files from Zenodo - https://doi.org/10.5281/zenodo.15161568. The key variables are: x0 (cross-shore position), z0 (initial bed level), z1 (final bed level), z0_rock (initial rock layer), z1_rock (final rock layer)
6. To view summary outputs for morphology classification in GIS, shapefiles are located in "data_1_shp.zip", for every 10'th profile (every 300 m alongshore).



FULL LIST OF ZIP FILES
- "_tutorial_scripts+data_22Jun2025.zip" --> Contains all Matlab (2024) scripts and cropped datasets required to run the profile preparation and classification algorithm and ShoreTrans model for a selected set of 7 profiles. For the tutorial scripts to run, the suite of ShoreTrans functions must be added to the Matlab path, these can be found in item "_ShoreTrans_functions_22Jun2025.zip" (listed below). The code can be run with a base Matlab installation, as long as the Profile Classificaiton and ShoreTrans functions are added to Matlab path).

- "_profile_prep+classify_functions_22Jun2025.zip" --> Matlab functions for Profile Preparation and Morphological Classification. These must be added to the Matlab path for the tutorial scripts to run.

- "_ShoreTrans_functions_22Jun2025.zip" --> ShoreTrans functions (update June 2025). These must be added to the Matlab path for the tutorial scripts to run.

- "data_1_shp.zip" --> VIC ShoreTrans model - Shapefile summary output of morphology classification.

- "data_2_run040_var_ST_good.zip" --> VIC ShoreTrans model - classification outputs and settings used in ShoreTrans for all profiles.

- "data_3_run040_var_z0.zip" --> VIC ShoreTrans model - input (initial, 2010) profiles (x0 = cross-shore distance; z0 = 2010 bed level; z0_rock = 2010 interpreted rock layer level).

- "data_4_run040_var_z1_hi.zip" --> VIC ShoreTrans model - output bed level (final, 2100) profiles (z1 = 2100 bed level; z1_lo = low bound uncertainty; z1_hi = high bound uncertainty).

- "data_5_run040_var_z1_hi.zip" --> VIC ShoreTrans model - output rock layer (final, 2100) profiles (z1_rock = 2100 rock level; z1_rock_lo = low bound uncertainty; z1_rock_hi = high bound uncertainty).

- "data_6_run040_STV.zip" --> VIC ShoreTrans model - short-term variability (storm erosion) for dune-backed profiles (output profiles: stv_z_lo - low end estimate; stv_z_av - mean estimate; stv_z_hi - high end estimate) 

- "working_code_prep_analysis.zip" --> Full set of scripts used for the Victorian ShoreTrans project and manuscript preparation. These are personal use scripts and other uses will encounter errors if they try to run them. Included to allow review.



















 





 














