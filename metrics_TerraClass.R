## Calculate metrics data for TerraBio for the sampling design using TerraClass
## data.

## This script is part of a series that calculates metrics for three separate
## maps, the 150 farm data from Imaflora (date unknown), MapBiomas 2018, and
## TerraClass 2014. The metrics are focused largely on forest cover as this is
## known to be a key landscape component for the butterflies and insects of
## interest in our eDNA sampling.

##
## Import --> [[Metrics]] --> Sampling
##

## ------- Load Packages -------------------------------------

library(raster)
library(landscapemetrics)
library(corrplot)
library("dplyr")
library("maptools")
library(landmetrics)

# Run this if you don't have the landmetrics package installed. library(remotes)
#
# remotes::install_github("jeffreyevans/landmetrics")
#
# The landmetrics package looks to be faster than the landscapemetrics package
# for doing focal window calculations. https://github.com/jeffreyevans/landmetrics


## ------- Load Data -----------------------------------------

source("import.R")

## Note, for landscapemetrics package, this must be as a raster with units in meters.

listlsm <- list_lsm()

## ------- Calculate Farm-based Metrics for TerraClass ----------------------------

## TerraClass data can be calculated both where farms are (using Imaflora
## boundaries) and using a moving window.

## Use the sample_lsm function to sample within farm boundaries.
## https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html

# 1 = Primary Natural Forest Vegetation
# 2 = Secondary Natural Forest Vegetation
# 3 = Silviculture
# 4 = Cultivated Pasture Cultivated Shrub
# 5 = Herbaceous Cultivated
# 6 = Perennial Agricultural Crop
# 7 = Seim-Perennial Agricultural Crop
# 8 = Temporary Agricultural Crop
# 9 = Mining
# 10 = Urbanized Area
# 11 = Other
# 12 = Not Observed
# 13 = Deforestation in the reference year
# 14 = Legal Amazon area that is not forest
# 15 = Bodies of water


## This is to recode 
isBecomes <- cbind(c(2),
                   c(1))

terraclass.2014.rast <- reclassify(terraclass.2014.rast, rcl=isBecomes)



## Edge complexity using Fractal Dimension Index (FDI), focus on forest patches.

FRACM.class.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "frac"#, 
        #all_classes = TRUE
    )

hist(FRACM.class.tc.30m$value[FRACM.class.tc.30m$class == 1 &
                                  FRACM.class.tc.30m$metric == "frac_mn"])


FRACM.landscape.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "frac"
    )

hist(FRACM.landscape.tc.30mm$value[FRACM.landscape.tc.30m$metric ==
                                       "frac_mn"])



# Total areas


TAREA.class.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "ca"
    )

hist(TAREA.class.tc.30m$value[TAREA.class.tc.30m$class == 1])


TAREA.landscape.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "ta"
    )

hist(TAREA.landscape.tc.30m$value, breaks = 30)




## Number of patches


NUMPAT.landscape.tc.30m <-
    sample_lsm(terraclass.2014.rast,
               y = farm.boundary.shp,
               level = "landscape",
               metric = "np")

hist(NUMPAT.landscape.tc.30m$value)



## Edge density

EDGE.class.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "ed"
    )

hist(EDGE.class.tc.30m$value[EDGE.class.tc.30m$class == 1])



EDGE.landscape.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "ed"
    )

hist(EDGE.landscape.tc.30m$value)



# patch density

PD.class.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "pd"
    )

hist(PD.class.tc.30m$value[PD.class.tc.30m$class == 1])


PD.landscape.tc.30m <-
    sample_lsm(terraclass.2014.rast,
               y = farm.boundary.shp,
               level = "landscape",
               metric = "pd")

hist(PD.landscape.tc.30m$value)


# percentage of landscape (PLAND)

PLAND.class.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "pland"
    )

hist(PLAND.class.tc.30m$value[PLAND.class.tc.30m$class == 1])




## Diversity of land covers--Shannon's diversity

SHDI.landscape.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "shdi"
    )

hist(SHDI.landscape.tc.30m$value)



## Contiguity Index


CONT.class.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "contig"
    )

hist(CONT.class.tc.30m$value[CONT.class.tc.30m$metric == "contig_mn" &
                                 CONT.class.tc.30m$class == 1])



CONT.landscape.tc.30m <-
    sample_lsm(
        terraclass.2014.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "contig"
    )

hist(CONT.landscape.tc.30m$value[CONT.landscape.tc.30m$metric == "contig_cv"])



## ------- Looking at farm-based correlation ---------------------------

## Method using corrplot. First, build the different class tables--not all farms have all classes, and asking landscapemetrics to include all patch types causes a vector size error 


# For the forest class
tc.forest <- tibble(
    plot_id = FRACM.class.tc.30m$plot_id[FRACM.class.tc.30m$class == 1 &
                                             FRACM.class.tc.30m$metric == "frac_mn"],
    FRACM.mean.forest = FRACM.class.tc.30m$value[FRACM.class.tc.30m$class == 1 &
                                                     FRACM.class.tc.30m$metric == "frac_mn"],
    TAREA.forest = TAREA.class.tc.30m$value[TAREA.class.tc.30m$class ==
                                                1],
    EDGE.forest = EDGE.class.tc.30m$value[EDGE.class.tc.30m$class == 1],
    PD.forest = PD.class.tc.30m$value[PD.class.tc.30m$class == 1],
    PLAND.forest = PLAND.class.tc.30m$value[PLAND.class.tc.30m$class == 1],
    CONT.forest = CONT.class.tc.30m$value[CONT.class.tc.30m$metric == "contig_mn" &
                                              CONT.class.tc.30m$class == 1]
    
)

corrplot::corrplot(cor(tc.forest[-1]), method = "shade")


tc.landscape <- tibble(
    plot_id = FRACM.landscape.tc.30m$plot_id[FRACM.landscape.tc.30m$metric=="frac_mn"],
    
    FRACM.mean.farm = FRACM.landscape.tc.30m$value[FRACM.landscape.tc.30m$metric=="frac_mn"],
    TAREA.farm = TAREA.landscape.tc.30m$value,
    NUMPAT.farm = NUMPAT.landscape.tc.30m$value,
    PD.farm = PD.landscape.tc.30m$value,
    EDGE.farm = EDGE.landscape.tc.30m$value,
    SHDI.farm = SHDI.landscape.tc.30m$value,
    CONT.mean.farm = CONT.landscape.tc.30m$value[CONT.landscape.tc.30m$metric == "contig_mn"],
    
)

corrplot(cor(tc.landscape[-1]), method = "shade")


# Combine

all.tc <- left_join(tc.landscape, tc.forest) 

all.tc[is.na(all.tc)] <- 0
corrplot(cor(all.tc[-1]), method = "shade")
pairs(all.tc[-1], lower.panel = panel.smooth, upper.panel = panel.cor)


## ------- Calculate moving window metrics for TerraClass ---------------------------
