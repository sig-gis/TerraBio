## This script takes eDNA data and environmental (mapping) data and combines it
## using species distribution models.

source("data_processing_IMPLEMENT.R")
library(plyr)
library(tidyr)
library(dplyr)

rm(list = ls()[grepl(ls(), pattern = "riaz|gill|zea")])

rare = 0


if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("sylvainschmitt/SSDM")


## ----- Habitat framework ------------------------------------------
# 
# Use hymenoptera & lepidoptera

# Need environmental variables--NDVI, tree height, and elevation along with land use should be good.

# Use SSDM Package
library(SSDM)

## Variable setup
envData <- load_var(path = "SmallEnvData/", tmp = TRUE, format = ".tif", files = "small_topoDiversity.tif")
envData

# Species data needs to be loaded with three columns: Species, Longitude, and Latitude. They should be in the same projection...

# For this, we need to have the longitude and latitude data and then match it up based on site.

motuPlots <-
    pivot_longer(
        invert_2_ReadOnly,
        cols = `SFX237-04P-03`:`SFX237-04P-02`,
        names_to = "plotName",
        values_to = "reads"
    )
motuPlots <- motuPlots[ motuPlots$reads > rare , ]
motuPlots <- motuPlots[ !grepl( "F-", motuPlots$plotName) , ] # remove Forest plots

specData <- left_join(motuPlots, sampledPlotsTable) %>% 
    left_join(keyInvertSpecies)
    
# Use the $order column so we can use hymenoptera and lepidoptera.
    
specData <- specData %>% dplyr::select(id, order, plotName, xCentroid, yCentroid) %>%
                  subset(specData$order %in% c('Lepidoptera', 'Hymenoptera'))
summary(specData)

## Models


# Used stacked with Lepidoptera and Hymenoptera

lepid_SDM <- modelling(
  algorithm = "SVM",
  Occurrences = specData[ specData$order == "Lepidoptera", ],
  Env = envData,
  Xcol = 'xCentroid',
  Ycol = 'yCentroid',
  verbose = TRUE,
  save = TRUE,
  ensemble.thresh = .1,
  cores = 8
)

knitr::kable(lepid_SDM@evaluation)
knitr::kable(lepid_SDM@variable.importance)
plot(lepid_SDM)
save(lepid_SDM, file = "Outputs/LEPID_SDM.RData")

# Used stacked with Lepidoptera and Hymenoptera

invert_SSDM <- stack_modelling(
  algorithms = "SVM",
  Occurrences = specData,
  Env = envData,
  Xcol = 'xCentroid',
  Ycol = 'yCentroid',
  verbose = TRUE,
  Spcol = "order",
  save = TRUE,
  ensemble.thresh = .1,
  cores = 8
)



invert_2_SSDM <-
  ensemble_modelling(
    c('RF'),
    subset(specData, specData$order %in% c('Lepidoptera', 'Hymenoptera')),
    Env = envData,
    Xcol = 'xCentroid',
    Ycol = 'yCentroid',
    verbose = TRUE,
    Spcol = "order",
    save = TRUE,
    ensemble.thresh = 0,
    cores = 8
  )
saveRDS(invert_2_SSDM, file = "invert_2_SSDM.RData")



plot(invert_2_SSDM@projection, main = 'SSDM\nfor Lepidoptera and Hymenoptera\nwith Maxent algorithm')


# evaluate the model
knitr::kable(invert_2_SSDM@evaluation)

# examine importance of environmental variables

knitr::kable(invert_2_SSDM@variable.importance)

  

# Possible packages:
# https://onlinelibrary.wiley.com/doi/10.1111/ecog.02671
# 
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.01881
# 
# https://rspatial.org/sdm/SDM.pdf
# 
#SSDM: https://mran.microsoft.com/web/packages/SSDM/vignettes/SSDM.html