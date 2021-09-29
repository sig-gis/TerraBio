## Calculate metrics data for TerraBio for the sampling design using MapBiomas
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


## ------- Load Data and Vars --------------------------------

source("import.R")

## Note, for landscapemetrics package, this must be as a raster with units in meters.

listlsm <- list_lsm()

cleanup <- TRUE             # True will remove intermediary processing steps for cleaner environment.
write.movingwindow <- FALSE  # True will cause moving window to run & write rasters. WARNING THIS IS VERY TIME INTENSIVE
mw.calcs <- FALSE # True will fix moving window to fill in spaces not calculated. It pulls from maps saved during write.movingwindow. 



## ------- Calculate Farm-based Metrics for MapBiomas ----------------------------

## MapBiomas data can be calculated both where farms are (using Imaflora
## boundaries) and using a moving window.

## Use the sample_lsm function to sample within farm boundaries.
## https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html

# 3 = Forest
# 4 = Savanna
# 5 = Mangrove
# 9 = Forest Plantation
# 10 = Non FOrest Natural FOrmation
# 11 = Wetland
# 12 = Grassland
# 13 = Other non forest natural formation
# 15 = Pasture
# 19 = Annual and Perennial crop
# 21 = Mosaic of Ag and Pasture
# 23 = Beach and Dune
# 24 = Urban INfrastructure
# 25 = other non-vegetated area
# 26 = water
# 27 = unknown
# 29 = rocky outcrop
# 30 = mining
# 32 = salt flat
# 33 = river lake ocean

## This is to recode if needed.
# isBecomes <- cbind(c(0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, etc, etc),
#                    c(0, 13, 12, 14, 14, 13, 14, 12, 14, 14, 13, 12))
# 
# mapbiomas.2018.rast <- reclassify(mapbiomas.2018.rast, rcl=isBecomes)



## Edge complexity using Fractal Dimension Index (FDI), focus on forest patches.

FRACM.class.mb.30m <-
    sample_lsm(
        mapbiomas.2018.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "frac"#, 
        #all_classes = TRUE
    )

hist(FRACM.class.mb.30m$value[FRACM.class.mb.30m$class == 3 &
                                    FRACM.class.mb.30m$metric == "frac_mn"])



# Total areas


TAREA.class.mb.30m <-
    sample_lsm(
        mapbiomas.2018.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "ca"
    )

hist(TAREA.class.mb.30m$value[TAREA.class.mb.30m$class == 3])


TAREA.class.mb.30m.buf <-
  sample_lsm(
    mapbiomas.2018.rast,
    y = farm.boundary.buffer,
    level = "class",
    metric = "ca"
  )

hist(TAREA.class.mb.30m.buf$value[TAREA.class.mb.30m.buf$class == 3])



## Number of patches


NUMPAT.class.mb.30m <-
    sample_lsm(mapbiomas.2018.rast,
               y = farm.boundary.shp,
               level = "class",
               metric = "np")

hist(NUMPAT.class.mb.30m$value[NUMPAT.class.mb.30m$class == 3])


NUMPAT.class.mb.30m.buf <-
  sample_lsm(mapbiomas.2018.rast,
             y = farm.boundary.buffer,
             level = "class",
             metric = "np")

hist(NUMPAT.class.mb.30m.buf$value[NUMPAT.class.mb.30m.buf$class == 3])


## Edge density

EDGE.class.mb.30m <-
    sample_lsm(
        mapbiomas.2018.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "ed"
    )

hist(EDGE.class.mb.30m$value[EDGE.class.mb.30m$class == 3])



# patch density

PD.class.mb.30m <-
    sample_lsm(
        mapbiomas.2018.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "pd"
    )

hist(PD.class.mb.30m$value[PD.class.mb.30m$class == 3])



# percentage of landscape

PLAND.class.mb.30m <-
    sample_lsm(
        mapbiomas.2018.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "pland"
    )

hist(PLAND.class.mb.30m$value[PLAND.class.mb.30m$class == 3])


PLAND.class.mb.30m.buf <-
  sample_lsm(
    mapbiomas.2018.rast,
    y = farm.boundary.buffer,
    level = "class",
    metric = "pland"
  )

hist(PLAND.class.mb.30m.buf$value[PLAND.class.mb.30m.buf$class == 3])



## Contiguity Index


CONT.class.mb.30m <-
    sample_lsm(
        mapbiomas.2018.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "contig"
    )

hist(CONT.class.mb.30m$value[CONT.class.mb.30m$metric == "contig_mn" &
                                   CONT.class.mb.30m$class == 3])




## Aggregation Index (AI)

AI.class.mb.30m <-
  sample_lsm(
    mapbiomas.2018.rast,
    y = farm.boundary.shp,
    level = "class",
    metric = "ai"
  )

hist(AI.class.mb.30m$value[ AI.class.mb.30m$class == 3])

AI.class.mb.30m.buf <-
  sample_lsm(
    mapbiomas.2018.rast,
    y = farm.boundary.buffer,
    level = "class",
    metric = "ai"
  )

hist(AI.class.mb.30m.buf$value[ AI.class.mb.30m.buf$class == 3])







## ------- Looking at farm-based correlation ---------------------------

## Method using corrplot. First, build the different class tables--not all farms have all classes, and asking landscapemetrics to include all patch types causes a vector size error 


# For the forest class
mb.forest <- tibble(
    plot_id = FRACM.class.mb.30m$plot_id[FRACM.class.mb.30m$class == 3 &
                                               FRACM.class.mb.30m$metric == "frac_mn"],
    FRACM.mean.forest = FRACM.class.mb.30m$value[FRACM.class.mb.30m$class == 3 &
                                                       FRACM.class.mb.30m$metric == "frac_mn"],
    TAREA.forest = TAREA.class.mb.30m$value[TAREA.class.mb.30m$class ==
                                                  3],
    NUMPAT.forest = NUMPAT.class.mb.30m$value[NUMPAT.class.mb.30m$class == 3],
    EDGE.forest = EDGE.class.mb.30m$value[EDGE.class.mb.30m$class == 3],
    PD.forest = PD.class.mb.30m$value[PD.class.mb.30m$class == 3],
    PLAND.forest = PLAND.class.mb.30m$value[PLAND.class.mb.30m$class == 3],
    CONT.forest = CONT.class.mb.30m$value[CONT.class.mb.30m$metric == "contig_mn" &
                                                CONT.class.mb.30m$class == 3],
    AI.forest = AI.class.mb.30m$value[ AI.class.mb.30m$class == 3]
    
)

cor(mb.forest)
corrplot::corrplot(cor(mb.forest[-1]), method = "shade")

# for the buffers
mb.buffer <- tibble(
    plot_id = TAREA.class.mb.30m.buf$plot_id[TAREA.class.mb.30m.buf$class == 3],
    
    TAREA.buffer = TAREA.class.mb.30m.buf$value[TAREA.class.mb.30m.buf$class == 3],
    NUMPAT.buffer = NUMPAT.class.mb.30m.buf$value[NUMPAT.class.mb.30m.buf$class == 3],
    PLAND.buffer = PLAND.class.mb.30m.buf$value[PLAND.class.mb.30m.buf$class == 3],
    AI.buffer = AI.class.mb.30m.buf$value[AI.class.mb.30m.buf$class == 3]
)

corrplot(cor(mb.buffer[-1]), method = "shade")


# Combine

all.mb <- left_join(mb.buffer, mb.forest)

all.mb[is.na(all.mb)] <- 0
corrplot(cor(all.mb[-1]), method = "shade")
pairs(all.mb[-1], lower.panel = panel.smooth, upper.panel = panel.cor)


## ------- Cleanup -------------------------------

if (cleanup == TRUE){

remove(FRACM.class.mb.30m, 
       TAREA.class.mb.30m, TAREA.class.mb.30m.buf,
       NUMPAT.class.mb.30m, NUMPAT.class.mb.30m.buf,
       EDGE.class.mb.30m,
       PD.class.mb.30m, 
       PLAND.class.mb.30m, PLAND.class.mb.30m.buf,
       CONT.class.mb.30m, 
       AI.class.mb.30m, AI.class.mb.30m.buf, 
       mb.forest, mb.buffer)

}




## ------- Calculate moving window metrics for MapBiomas ---------------------------

## implementation with landscape metrics

## only landscape metrics are allowed for moving windows.

## Note--this is prohibitively time intensive currently. 

# moving_window <- matrix(1, nrow = 25, ncol = 25)
# 
# landscape.metrics <-
#     c("lsm_l_frac_mn",
#       "lsm_l_np",
#       "lsm_l_ed",
#       "lsm_l_pd",
#       "lsm_l_shdi",
#       "lsm_l_contig")
# 
# window.landscape.mb <- window_lsm(mapbiomas.2018.rast, 
#                               window = moving_window,
#                               what = "lsm_l_pd", progress = T
#                               )


## Implementation with landmetrics: metrics available:
## https://github.com/jeffreyevans/landmetrics/blob/master/R/focal.lmetrics.R
## This package allows for class type calculations

# Here are the functions: note that Shannon's and Contiguity aren't available.
# mean.frac.dim.index - mean of fractal dimension index.
# mean.patch.area - average area of patches.
# total.area - the sum of the areas (m2) of all patches of the corresponding patch type.
# n.patches - the number of patches of a particular patch type or in a class.
# edge.density - edge length on a per unit area basis that facilitates comparison among landscapes of varying size.
# patch.density - the numbers of patches of the corresponding patch type divided by total landscape area (m2).
# prop.landscape - the proportion of the total landscape represented by this class
# patch.cohesion.index - measures the physical connectedness of the corresponding patch type.
# aggregation.index - computed simply as an area-weighted mean class aggregation index, where each class is weighted by its proportional area in the landscape.
# effective.mesh.size - equals 1 divided by the total landscape area (m2) multiplied by the sum of patch area (m2) squared, summed across all patches in the landscape.

## The if condition here means that these processing intensive calculations
## occur ONLY when needed.

if (write.movingwindow == TRUE){

FRAC.forest.window.mb <-
    focal.lmetrics(
        mapbiomas.2018.rast,
        w = 5,
        land.value = 3,
        metric = "mean.frac.dim.index"
    )

writeRaster(FRAC.forest.window.mb, "../Xingue basin/MapBiomas/FRAC.forest.5window.mb.tif")

TAREA.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 5,
    land.value = 3,
    metric = "total.area"
  )

writeRaster(TAREA.forest.window.mb, "../Xingue basin/MapBiomas/TAREA.forest.5window.mb.tif")


NUMPAT.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 5,
    land.value = 3,
    metric = "n.patches"
  )

writeRaster(NUMPAT.forest.window.mb, "../Xingue basin/MapBiomas/NUMPAT.forest.5window.mb.tif")

EDGE.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 5,
    land.value = 3,
    metric = "edge.density"
  )

writeRaster(EDGE.forest.window.mb, "../Xingue basin/MapBiomas/EDGE.forest.5window.mb.tif")
 
PD.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 5,
    land.value = 3,
    metric = "patch.density"
  )

writeRaster(PD.forest.window.mb, "../Xingue basin/MapBiomas/PD.forest.5window.mb.tif")


PLAND.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 5,
    land.value = 3,
    metric = "prop.landscape"
  )

writeRaster(PLAND.forest.window.mb, "../Xingue basin/MapBiomas/PLAND.forest.5window.mb.tif")

PLAND.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 15,
    land.value = 3,
    metric = "prop.landscape"
  )

writeRaster(PLAND.forest.window.mb, "../Xingue basin/MapBiomas/PLAND.forest.15window.mb.tif")

COHESION.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 5,
    land.value = 3,
    metric = "patch.cohesion.index"
  )

writeRaster(COHESION.forest.window.mb, "../Xingue basin/MapBiomas/COHESION.forest.5window.mb.tif")


AI.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 5,
    land.value = 3,
    metric = "aggregation.index"
  )

writeRaster(AI.forest.window.mb, "../Xingue basin/MapBiomas/AI.forest.5window.mb.tif")

AI.forest.window.mb <-
  focal.lmetrics(
    mapbiomas.2018.rast,
    w = 15,
    land.value = 3,
    metric = "aggregation.index"
  )

writeRaster(AI.forest.window.mb, "../Xingue basin/MapBiomas/AI.forest.15window.mb.tif")


}


## -------- Modify moving window matrix ---------------------------

## The moving window rasters that get computed are only computed at the edge of
## forests. So we need to fill in the na values.


if(mw_calcs == TRUE){

## __FRAC__

## load the data
  FRAC.forest.window.mb <- raster("../Xingue basin/MapBiomas/FRAC.forest.5window.mb.tif")

## assign desired values for forest and non-forest
  forest = 1
  nonforest = 1

## Get rid of NA in the data
  FRAC.forest.window.mb[is.na(FRAC.forest.window.mb)] <- 999
  
## create the base map based on the classified raster
  temp.forest <- mapbiomas.2018.rast
  temp.forest[temp.forest != 3] <- forest
  temp.forest[temp.forest == 3] <- nonforest


## Sum the two maps
  
  FRAC.forest.window.mb <-
    overlay(
      FRAC.forest.window.mb,
      temp.forest,
      fun = function(x, y) {
        ifelse(x == 999 & y == nonforest, nonforest, ifelse(x == 999 & y == forest, forest, x))
      }
    )
  
  writeRaster(
    FRAC.forest.window.mb,
    "../Xingue basin/MapBiomas/FRAC.forest.5window.mb.filled.tif"
  )
      

## __TAREA__
  
## load the data
  TAREA.forest.window.mb <- raster("../Xingue basin/MapBiomas/TAREA.forest.5window.mb.tif")
  
## assign desired values for forest and non-forest
  forest = 22500
  nonforest = 0
  
## Get rid of NA in the data
  TAREA.forest.window.mb[is.na(TAREA.forest.window.mb)] <- 999
  
## create the base map based on the classified raster
  temp.forest <- mapbiomas.2018.rast
  temp.forest[temp.forest != 3] <- forest
  temp.forest[temp.forest == 3] <- nonforest
  
  
## Sum the two maps
  
  TAREA.forest.window.mb <-
    overlay(
      TAREA.forest.window.mb,
      temp.forest,
      fun = function(x, y) {
        ifelse(x == 999 & y == nonforest, nonforest, ifelse(x == 999 & y == forest, forest, x))
      }
    )
  
  writeRaster(
    TAREA.forest.window.mb,
    "../Xingue basin/MapBiomas/TAREA.forest.5window.mb.filled.tif",
    overwrite = TRUE
  )
  
## __EDGE__
  
## load the data
  EDGE.forest.window.mb <- raster("../Xingue basin/MapBiomas/EDGE.forest.5window.mb.tif")
  
## assign desired values for forest and non-forest
  forest = 0
  nonforest = 0
  
## Get rid of NA in the data
  EDGE.forest.window.mb[is.na(EDGE.forest.window.mb)] <- 999
  
## create the base map based on the classified raster
  temp.forest <- mapbiomas.2018.rast
  temp.forest[temp.forest != 3] <- forest
  temp.forest[temp.forest == 3] <- nonforest
  
  
## Sum the two maps
  
  EDGE.forest.window.mb <-
    overlay(
      EDGE.forest.window.mb,
      temp.forest,
      fun = function(x, y) {
        ifelse(x == 999 & y == nonforest, nonforest, ifelse(x == 999 & y == forest, forest, x))
      }
    )
  
  writeRaster(
    EDGE.forest.window.mb,
    "../Xingue basin/MapBiomas/EDGE.forest.5window.mb.filled.tif"
  )
  
## __PLAND__
  
## load the data
  PLAND.forest.window.mb <- raster("../Xingue basin/MapBiomas/PLAND.forest.5window.mb.tif")
  
## assign desired values for forest and non-forest
  forest = 1
  nonforest = 0
  
## Get rid of NA in the data
  PLAND.forest.window.mb[is.na(PLAND.forest.window.mb)] <- 999
  
## create the base map based on the classified raster
  temp.forest <- mapbiomas.2018.rast
  temp.forest[temp.forest != 3] <- forest
  temp.forest[temp.forest == 3] <- nonforest
  
  
## Sum the two maps
  
  PLAND.forest.window.mb <-
    overlay(
      PLAND.forest.window.mb,
      temp.forest,
      fun = function(x, y) {
        ifelse(x == 999 & y == nonforest, nonforest, ifelse(x == 999 & y == forest, forest, x))
      }
    )
  
  writeRaster(
    PLAND.forest.window.mb,
    "../Xingue basin/MapBiomas/PLAND.forest.5window.mb.filled.tif"
  )
  
## __COHESION__
  
  ## load the data
  COHESION.forest.window.mb <- raster("../Xingue basin/MapBiomas/COHESION.forest.5window.mb.tif")
  
  ## assign desired values for forest and non-forest
  forest = 8
  nonforest = 0
  
  ## Get rid of NA in the data
  COHESION.forest.window.mb[is.na(COHESION.forest.window.mb)] <- 999
  
  ## create the base map based on the classified raster
  temp.forest <- mapbiomas.2018.rast
  temp.forest[temp.forest != 3] <- forest
  temp.forest[temp.forest == 3] <- nonforest
  
  
  ## Sum the two maps
  
  COHESION.forest.window.mb <-
    overlay(
      COHESION.forest.window.mb,
      temp.forest,
      fun = function(x, y) {
        ifelse(x == 999 & y == nonforest, nonforest, ifelse(x == 999 & y == forest, forest, x))
      }
    )
  
  writeRaster(
    COHESION.forest.window.mb,
    "../Xingue basin/MapBiomas/COHESION.forest.5window.mb.filled.tif"
  )
  
  
## __AI__
  
  ## load the data
  AI.forest.window.mb <- raster("../Xingue basin/MapBiomas/AI.forest.5window.mb.tif")
  
  ## assign desired values for forest and non-forest
  forest = 100
  nonforest = 0
  
  ## Get rid of NA in the data
  AI.forest.window.mb[is.na(AI.forest.window.mb)] <- 999
  
  ## create the base map based on the classified raster
  temp.forest <- mapbiomas.2018.rast
  temp.forest[temp.forest != 3] <- forest
  temp.forest[temp.forest == 3] <- nonforest
  
  
  ## Sum the two maps
  
  AI.forest.window.mb <-
    overlay(
      AI.forest.window.mb,
      temp.forest,
      fun = function(x, y) {
        ifelse(x == 999 & y == nonforest, nonforest, ifelse(x == 999 & y == forest, forest, x))
      }
    )
  
  writeRaster(
    AI.forest.window.mb,
    "../Xingue basin/MapBiomas/AI.forest.5window.mb.filled.tif"
  )
  
}