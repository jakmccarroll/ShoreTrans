# ShoreTrans (June 2025 update) **
This update contains ShoreTrans functions and Profile Morphological Classification functions associated with the article:

McCarroll, R. J., Kennedy, D. M., & Ierodiaconou, D. (2025). Morphologically adaptive modelling of sea level rise induced coastal erosion impacts for south-East Australia. Marine Geology, 107602.

DOI: https://doi.org/10.1016/j.margeo.2025.107602

Email: jak.mccarroll@gmail.com 

ZENODO REPOSITORY

The ShoreTrans and Morphological Classification functions as well as input / output datasets are available in a Zenodo repository: https://doi.org/10.5281/zenodo.15161568


QUICK INSTRUCTIONS

1. Download and unzip "_tutorial_scripts+data_22Jun2025.zip"
2. Download and unzip "_ShoreTrans_functions_22Jun2025.zip" ---> add all functions to the Matlab path
3. Download and unzip "_profile_prep+classify_functions_22Jun2025.zip" ---> add all functions to the Matlab path
3. Open script "../tutorial/scripts/tut01_prep_classify.m" --> read through instructions, run the profile classification for the 7 example profiles
4. Open script "../tutorial/scripts/tut02_ST_run.m" --> run ShoreTrans for the 7 example profiles


+++++++++++

TITLE: Morphologically adaptive modelling of sea level rise induced coastal erosion impact for south-east Australia

AUTHORS: Robert Jak McCarroll1, David M. Kennedy1, Daniel Ierodiaconou2

1 School of Geography, Earth and Atmospheric Sciences, The University of Melbourne, Parkville, Victoria, Australia

2 School of Life and Environmental Sciences, Deakin University, Warrnambool, Victoria, Australia

ABSTRACT

Sea level rise induced coastal erosion represents an impending threat to the world's coastlines. A critical control on coastal recession is the onshore accommodation space available for a receding beach to occupy. Despite the importance of morphologic controls, coastal change models applied at regional-scale over ~100-year time-frames typically address a limited range of coastal morphologies.

ShoreTrans is a shoreface translation model that kinematically projects profile change, based on inputs of sea level rise and sediment budget imbalances. This work presents updates to enable broad-scale model application (1000's km), automating classification and adapting to the morphology of individual profiles, including: (1) dunes; (2) cliffs; (3) bluffs and ridges; (4) inter-subtidal rock outcrops; (5) protection structures; (6) low-energy environments; and (7) short-term dune erosion. The model was applied to Victoria, Australia (2000 km coastline, 30 m spaced transects), for a scenario of 1 m sea level rise from 2010 to 2100, using a single time step and simplified treatment of uncertainty and sediment budget. Shoreline trends and variability were determined from satellite extracted shorelines.

Mean projected shoreline recession of 43 m is 20 % lower than for a simple parameterization (uncertainty range 58 % lower to 8 % higher), due to sediment transfer from the backshore to the active shoreface (e.g., dune encroachment) and hard backshores restricting shoreline movements (e.g., cliffs, seawalls). Low dunes exhibited the highest recession rates, due to rollover (68 m long-term recession). Total setback extent, including short-term variability, is projected to exceed 182 m in 5 % of low dune areas. High rates of beach loss were associated with beaches fronting hard cliffs (55 % beach loss) and seawalls (80 %). The worst impacts are expected for rocky, sediment poor coastlines, such as the Great Ocean Road Surf Coast, where a loss of 30 % to 50 % of beaches is projected, not accounting for infrastructure and potential management interventions.

Automated morphological adaptation represents a step-change for regional scale coastal change assessment. The method also allows for coupling of future erosion to inundation hazards, by interpolating a 3D surface of future morphology. At local-scale, ShoreTrans is suitable to add to hybrid models, providing a means to improve future coastal change projections.

22/6/2025

+++++++++++
# ShoreTrans (2021 version)

ShoreTrans: a simple, rules-based, user-input driven, shoreface translation and sediment budgeting model, that applies the surveyed 2D-profile (not a parameterization), for estimating change to realistic coastlines, resulting from sea level rise and variations in sediment supply, while accounting for armouring, hard-rock cliffs and outcropping rocks.  The tool can be applied to sand, gravel, rock and engineered coasts at a temporal scale of 10–100 years, accounting for shoreline trends as well as variability. The method accounts for: (1) dune encroachment/accretion; (2) barrier rollback; (3) non-erodible layers; (4) seawalls; (5) lower shoreface transport; (6) alongshore rotation; and (7) other sources and sinks. 

The ShoreTrans model is associated with the article:
"A rules-based shoreface translation and sediment budgeting tool for estimating coastal change: ShoreTrans" McCarroll et al., 2021, Marine Geology.
https://www.sciencedirect.com/science/article/pii/S0025322721000487

Code is in Matlab.

4/5/2021
