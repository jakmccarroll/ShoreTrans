**GITHUB - SHORETRANS (June 2025) **
This update contains ShoreTrans functions and Profile Morphological Classification functions associated with the article:
McCarroll, R. J., Kennedy, D. M., & Ierodiaconou, D. (2025). Morphologically adaptive modelling of sea level rise induced coastal erosion impacts for south-East Australia. Marine Geology, 107602.

DOI: https://doi.org/10.1016/j.margeo.2025.107602

Email: jak.mccarroll@gmail.com 

ZENODO REPOSITORY
The ShoreTrans and Morphological Classification functions as well as input / output datasets are available in a Zenodo repository
ZENODO LINK: https://doi.org/10.5281/zenodo.15161568


QUICK INSTRUCTIONS
1. Download and unzip "_tutorial_scripts+data_22Jun2025.zip"
2. Download and unzip "_ShoreTrans_functions_22Jun2025.zip" ---> add all functions to the Matlab path
3. Download and unzip "_profile_prep+classify_functions_22Jun2025.zip" ---> add all functions to the Matlab path
3. Open script "../tutorial/scripts/tut01_prep_classify.m" --> read through instructions, run the profile classification for the 7 example profiles
4. Open script "../tutorial/scripts/tut02_ST_run.m" --> run ShoreTrans for the 7 example profiles


+++++++++++
# ShoreTrans (2021 version)

ShoreTrans: a simple, rules-based, user-input driven, shoreface translation and sediment budgeting model, that applies the surveyed 2D-profile (not a parameterization), for estimating change to realistic coastlines, resulting from sea level rise and variations in sediment supply, while accounting for armouring, hard-rock cliffs and outcropping rocks.  The tool can be applied to sand, gravel, rock and engineered coasts at a temporal scale of 10–100 years, accounting for shoreline trends as well as variability. The method accounts for: (1) dune encroachment/accretion; (2) barrier rollback; (3) non-erodible layers; (4) seawalls; (5) lower shoreface transport; (6) alongshore rotation; and (7) other sources and sinks. 

The ShoreTrans model is associated with the article:
"A rules-based shoreface translation and sediment budgeting tool for estimating coastal change: ShoreTrans" McCarroll et al., 2021, Marine Geology.
https://www.sciencedirect.com/science/article/pii/S0025322721000487

Code is in Matlab.

Some elements presented in the article are still to be added or need testing/improvement, including:
- Time series of shoreface translation with SLR curve
- Probabilistic sampling and Monte Carlo simulations
- Non-erodible layer (rock outcrops) needs further testing
- Optimization routine needs re-writing

jak.mccarroll@plymouth.ac.uk

Jak McCarroll, 4/5/2021
