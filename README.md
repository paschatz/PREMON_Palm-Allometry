[![DOI](https://zenodo.org/badge/600030070.svg)](https://zenodo.org/doi/10.5281/zenodo.10390497)

# Height-diameter allometry for a dominant palm to improve understanding of carbon and forest dynamics in forests of Puerto Rico

This is a publicly accessible implementation of our study, providing the R-script containing all statistical analyses, figures and tables. To run the code in a reproducible envoroment we recommend to run code through ```PREMON_allometry.Rproj```. This will set the working directory to the root of the project.

## Abstract:
Tropical forests are major components of Earth’s carbon stocks, but their diversity and structural complexity challenge our ability to accurately estimate carbon stocks and dynamics. Palms, in particular, are prominent components of many tropical forests that have unique anatomical, physiological, and allometric differences from dicot trees, which impede accurate estimates of their above-ground biomass (AGB). In this study, we focused on improving height estimates and, ultimately, AGB estimates for a highly abundant palm in Puerto Rico, *Prestoea acuminata*. Based on field measurements of 1,215 individuals, we found a strong relationship between stem height and diameter at breast height. We also found evidence that height-diameter allometry of *P. acuminata* is mediated by various sources of environmental heterogeneity such as elevation, slope, and neighborhood crowding. We then examined variability in palm AGB estimates derived from two different models developed to estimate palm AGB. Finally, we applied our novel H:D allometric model with a site- and species-specific AGB model to hindcast palm AGB dynamics in the Luquillo Forest Dynamics Plot during a 35-year period of post-hurricane recovery. Overall, our study provides improved estimates of AGB in wet forests of Puerto Rico and will facilitate novel insights to the dynamics of palms in tropical forests.

## Keywords: 
Aboveground biomass, Arecaceae, environmental heterogeneity, forest dynamics, neotropics, *Prestoea acuminata var. montana* (Graham)

## Citations:

To cite this repository:

- Chatzopoulos, P., & Muscarella, R. (2023). paschatz/PREMON_Palm-Allometry: Height-diameter allometry for a dominant palm to improve understanding of carbon and forest dynamics in forests of Puerto Rico (v2.0.0). Zenodo. https://doi.org/10.5281/zenodo.10417497

To cite our data:

- Muscarella, R., P. Chatzopoulos, and R. Lammerant. 2023. Measurements of height, diameter at breast height and basal diameter for *Prestoea acuminata* at the Luquillo Forest Dynamics Plot (LFDP), Puerto Rico in January 2020. ver 2. Environmental Data Initiative. https://doi.org/10.6073/pasta/7c635efb584b1d254d8b713361b13ad0 (Accessed 2023-12-21).

### Folder description:
- ```clean_data``` = Here you can find the clean data that we used for our analysis. This data are an export of ```scripts/PREMON_data_prep.r``` script.

- ```figures``` = Here you can find the exported figures from the ```scripts/PREMON_analysis_script.R``` script.

- ```scripts``` = Here you can find the R scripts that we used for our analysis.

- ```tables``` = Here you can find the exported tables from the ```scripts/PREMON_analysis_script.R``` script.

- ```raw_data``` = Here you can find the raw data that we used for our analysis.

- ```renv``` = Here is stored a reproducible R environemnt. Users can reproduce the exact versions of r-packages we used during the analysis. This is an export of [renv](https://rstudio.github.io/renv/articles/renv.html) package.
