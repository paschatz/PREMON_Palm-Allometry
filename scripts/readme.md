```scripts/PREMON_data_prep.R```contains the script for cleaning and merging raw data files to prepare them for analysis. 
We used the EDIutils, readr, and tidyverse libraries to read and manipulate data. 
```"palm_allometry.csv,"``` contains the raw data with the measurements of *P. acuminata* individuals from January 2020.
To adress Q2 of the manuscript, we joined our dataset with the LFDP census data (census 6) and LFDP physical environment data from specific EDI repositories. 
We exported a final dataset where we stored our measured palms with the physical environment attributes for each stem,
such as elevation, slope, and the crowding index. Additionally, we calculated the slenderness ratio (SR) for each palm. 
SR was used as response variable for the multivariete regration in Q2. Finally, we filtered out individuals with height less
than 1.3 m or dbh measured at a different height (than 1.30 m). These are also the conditions under which we parameterized the H:D model in Q1. 
The resulting cleaned dataset is then saved as "raw_data/Data_for_analysis.csv" and used in ```PREMON_analysis.R```.

```PREMON_analysis.R``` contains all the code with which we generated the results of our manuscript. It is divided into sections that follow the specific four research questions of the manuscript. 

# Packages used:

- Attali D, Baker C (2023). _ggExtra: Add Marginal Histograms to
  'ggplot2', and More 'ggplot2' Enhancements_. R package version 0.10.1,
  <https://CRAN.R-project.org/package=ggExtra>.
  
- Kassambara A (2023). _ggpubr: 'ggplot2' Based Publication Ready
  Plots_. R package version 0.6.0,
  <https://CRAN.R-project.org/package=ggpubr>.

- Ogle DH, Doll JC, Wheeler AP, Dinno A (2023). _FSA: Simple Fisheries
  Stock Assessment Methods_. R package version 0.9.5,
  <https://CRAN.R-project.org/package=FSA>.

- R Core Team (2023). _R: A Language and Environment for Statistical
  Computing_. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.

- Rejou-Mechain M, Tanguy A, Piponiot C, Chave J, Herault B (2017).
  “BIOMASS : an R package for estimating above-ground biomass and its
  uncertainty in tropical forests.” _Methods in Ecology and Evolution_,
  *8*(9). ISSN 2041210X, doi:10.1111/2041-210X.12753
  <https://doi.org/10.1111/2041-210X.12753>.

- Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco
  Sciaini, and Cédric Scherer (2023). viridis(Lite) -
  Colorblind-Friendly Color Maps for R. viridis package version 0.6.4.

- Smith C (2023). _EDIutils: An API Client for the Environmental Data
  Initiative Repository in R_. <https://github.com/ropensci/EDIutils>.

- Ushey K, Wickham H (2023). _renv: Project Environments_. R package
  version 1.0.3, <https://CRAN.R-project.org/package=renv>.

- Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R,
  Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller
  E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V,
  Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to
  the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686.
  doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.





  











