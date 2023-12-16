####  Clean and merge raw data files for analysis data
library(EDIutils)
library(readr)
library(tidyverse)

## Read field data
data_raw <- read_csv("raw_data/palm_allometry.csv")

# Data for manipulation
data <- data_raw

# Remove unused columns
data$second.coord <- NULL
data$exposure <- NULL
data$h_dbh[data$height >= 1.3 & is.na(data$h_dbh)] <- 1.30

# "palm_allometry.csv" contains measurements of <i>P.acuminata</i> from January 2019. 
# Tag = Unique tag number for the individual stem
# basal_d = diameter at base of stem (cm)
# dbh = diameter at breast height (cm) 
# height = stem height (m)
# h_dbh = height at which dbh was measured (m)

## Read LFDP census data (census 6) to get quadrat of each (tagged) palm
# Original data and packageId information here:
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-luq.119.1545979

# Access raw data
raw6 <- EDIutils::read_data_entity(packageId = "knb-lter-luq.119.1545979", 
                                   entityId = "325c43057e0dd4e1cd6a13fa5125a76d")

# Save raw data in machine readable format
census6 <- readr::read_csv(file = raw6)
data$quadrat <- census6$Quadrat[match(data$Tag, census6$Tag)]

## Read LFDP physical environment data
# Original data and packageId information here:
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-luq.47.381052

rawenv <- EDIutils::read_data_entity(packageId = "knb-lter-luq.47.381050",
                                     entityId = "21c7cce729a96de675aaa9281526eb4a")
env <- readr::read_csv(file = rawenv)

## Read LFDP neighborhood crowding index data (computed separately from full LFDP data for census 6)
nci <- readRDS("raw_data/nci_data.RDA")
nci <- nci[nci$Census == 6,]

## Bind field data with physical environment dataset
data$Elev <- env$Elev[match(data$quadrat, env$Quadrat)]
data$Slope <- env$Slope[match(data$quadrat, env$Quadrat)]
data$nci <- nci$nci[match(data$Tag, nci$Tag)]

# Calculate slenderness ratio (SR)
data$SR <- data$height/(data$dbh*0.01)

# write.csv(data, "Data_for_analysis-including-small-palms.csv", row.names = F)

# Remove palms shorter than 1.3 m and where dbh was measured somewhere else than 1.3 m
data <- data[data$height >= 1.3 &  data$h_dbh == 1.30,]

write.csv(data, "clean_data//Data_for_analysis.csv", row.names = FALSE)