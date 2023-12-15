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
