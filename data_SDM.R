## This script takes eDNA data and environmental (mapping) data and combines it
## using species distribution models.


source("data_processing_IMPLEMENT.R")
library(tidyr)
library(dplyr)
#install.packages(c('raster', 'rgdal', 'dismo', 'rJava'))
library(maptools)
library(dismo)
library(raster)
library(rJava)
library(rgdal)

rm(list = ls()[grepl(ls(), pattern = "riaz|gill|zea")])

rare = 0
SDM_species <- c("Hymenoptera", "Lepidoptera")


## GIS data for the landscape

TCH2019 <- raster("EnvironmentalData/tch2019.tif")
topoDiversity <- raster("EnvironmentalData/topoDiversity.tif")
ndvi <- raster("EnvironmentalData/ndvi_4326.tif")
classification <- raster("EnvironmentalData/classification_20m.tif")


## ---- Subset Species Data -----------------------------------------

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
  left_join(keyInvertSpecies) %>% 
  dplyr::select(id, order, plotName, xCentroid, yCentroid)

# Use hymenoptera & lepidoptera
specData <- specData[ specData$order %in% SDM_species, ]

specData_LL <- dplyr::select(specData, xCentroid, yCentroid)


# Plot
# Determine geographic extent of our data
max.lat <- ceiling(max(specData$yCentroid))
min.lat <- floor(min(specData$yCentroid))
max.lon <- ceiling(max(specData$xCentroid))
min.lon <- floor(min(specData$xCentroid))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

# Load the data to use for our base map
data(wrld_simpl)

# Plot the base map
plot(wrld_simpl, 
     xlim = c(min.lon-1, max.lon+1),
     ylim = c(min.lat-1, max.lat+1),
     axes = TRUE, 
     col = "grey95")

# Add the points for individual observation
points(x = specData$xCentroid, 
       y = specData$yCentroid, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)
# And draw a little box around the graph
box()


## ----- Environmental data ---------------------------------------------
# Need environmental variables--NDVI, tree height, and elevation along with land use should be good.
## Variable setup
classificationCrop <- crop(x = classification, y = geographic.extent)


# test model
test.model <- bioclim(x = ndvi, p = specData_LL)

test.predict <- dismo::predict(object = test.model,
               x = ndvi,
               ext = geographic.extent)


## Models

# Used stacked with Lepidoptera and Hymenoptera
# Use maxent because it does as well as ensemble.

invert_2_SSDM <- modelling('MAXENT', subset(specData, specData$order %in% c('Lepidoptera', 'Hymenoptera')),
                          envData, Xcol = 'xCentroid', Ycol = 'yCentroid', verbose = FALSE)
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