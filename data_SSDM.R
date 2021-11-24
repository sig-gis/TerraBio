## This script takes eDNA data and environmental (mapping) data and combines it
## using species distribution models.

source("data_processing_IMPLEMENT.R")
library(plyr)
library(tidyr)
library(dplyr)

rm(list = ls()[grepl(ls(), pattern = "riaz|gill|zea")])

rare = 0


# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# 
# devtools::install_github("sylvainschmitt/SSDM")


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

## ----- Models ---------------------------------------------------------
#simple
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


# more environmental data

envData <- load_var(path = "SmallEnvData/", tmp = TRUE, format = ".tif", categorical = "small_classification_20m")
envData
plot(envData)

lepid_4ed_SDM <- modelling(
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

knitr::kable(lepid_4ed_SDM@evaluation)
knitr::kable(lepid_4ed_SDM@variable.importance)
plot(lepid_4ed_SDM)
save(lepid_4ed_SDM, file = "Outputs/LEPID_4ed_SDM.RData")

# try ensemble
lepid_4ed_ensemb <- ensemble_modelling(
  algorithms = c("SVM", "GLM", "RF"),
  Occurrences = specData[ specData$order == "Lepidoptera", ],
  Env = envData,
  Xcol = 'xCentroid',
  Ycol = 'yCentroid',
  verbose = TRUE,
  save = TRUE,
  ensemble.thresh = .5,
  cores = 8, name = "Lepidoptera"
)

# hymenoptera ensemble
hymen_4ed_ensemb <- ensemble_modelling(
  algorithms = c("SVM", "GLM", "RF"),
  Occurrences = specData[ specData$order == "Hymenoptera", ],
  Env = envData,
  Xcol = 'xCentroid',
  Ycol = 'yCentroid',
  verbose = TRUE,
  save = TRUE,
  ensemble.thresh = .5,
  cores = 8, name = "Hymenoptera"
)

# CTA,MARS gives poor results, as do GBM, ANN, 
#GAM does not work

SSDM_4ed_ensemb <- stacking(lepid_4ed_ensemb, hymen_4ed_ensemb)

save.stack(SSDM_4ed_ensemb)

plot(SSDM_4ed_ensemb, main = 'SSDM\nfor Lepidoptera and Hymenoptera\nwith Maxent algorithm')

# # evaluate the model
# knitr::kable(invert_2_SSDM@evaluation)
# 
# # examine importance of environmental variables
# 
# knitr::kable(invert_2_SSDM@variable.importance)

## ----- try full map --------------------------------------
## Variable setup
envDataBig <- load_var(path = "EnvironmentalData/", tmp = TRUE, format = ".tif", categorical = "classification_20m")
envDataBig

## Lepidoptera
lepid_big_ensemb <- ensemble_modelling(
  algorithms = c("SVM", "GLM"),
  Occurrences = specData[ specData$order == "Lepidoptera", ],
  Env = envDataBig,
  Xcol = 'xCentroid',
  Ycol = 'yCentroid',
  verbose = TRUE,
  save = TRUE,
  ensemble.thresh = .5,
  cores = 8, name = "LepidopteraBig"
)

## Hymenoptera
lepid_big_ensemb <- ensemble_modelling(
  algorithms = c("SVM", "GLM"),
  Occurrences = specData[ specData$order == "Hymenoptera", ],
  Env = envDataBig,
  Xcol = 'xCentroid',
  Ycol = 'yCentroid',
  verbose = TRUE,
  save = TRUE,
  ensemble.thresh = .5,
  cores = 8, name = "HymenopteraBig"
)




# Possible packages:
# https://onlinelibrary.wiley.com/doi/10.1111/ecog.02671
# 
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.01881
# 
# https://rspatial.org/sdm/SDM.pdf
# 
#SSDM: https://mran.microsoft.com/web/packages/SSDM/vignettes/SSDM.html