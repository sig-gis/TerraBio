## Calculate metrics data for TerraBio for the sampling design using Thibaud's
## data.

## This script is part of a series that calculates metrics for multiple maps.
## The metrics are focused largely on forest cover as this is known to be a key
## landscape component for the butterflies and insects of interest in our eDNA
## sampling.

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
write.movingwindow <- TRUE  # True will cause moving window to run & write rasters. WARNING THIS IS VERY TIME INTENSIVE




## ------- Calculate Farm-based Metrics for Thibaud ----------------------------

## Thibaud's data can be calculated both where farms are (using Imaflora
## boundaries) and using a moving window.

## Use the sample_lsm function to sample within farm boundaries.
## https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html

# -1 = No Data
# 0 = Savanna/bushes/scattered trees
# 1 = Grassland
# 2 = Mine
# 3 = Native
# 4 = Season
# 5 = Urban
# 6 = Water
# 7 = Cocoa


## This is to recode if needed.
# isBecomes <- cbind(c(0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, etc, etc),
#                    c(0, 13, 12, 14, 14, 13, 14, 12, 14, 14, 13, 12))
# 
# thibaud.2020.rast <- reclassify(thibaud.2020.rast, rcl=isBecomes)


## This is all commented out b/c Thibaud's layer doesn't match up with the
## farms--that one to the west is the problem.

# 
# ## Edge complexity using Fractal Dimension Index (FDI), focus on forest patches.
# 
# FRACM.class.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "frac"#, 
#         #all_classes = TRUE
#     )
# 
# hist(FRACM.class.thib.30m$value[FRACM.class.thib.30m$class == 3 &
#                                   FRACM.class.thib.30m$metric == "frac_mn"])
# 
# 
# FRACM.landscape.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "landscape",
#         metric = "frac"
#     )
# 
# hist(FRACM.landscape.thib.30m$value[FRACM.landscape.thib.30m$metric ==
#                                       "frac_mn"])
# 
# 
# 
# # Total areas
# 
# 
# TAREA.class.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "ca"
#     )
# 
# hist(TAREA.class.thib.30m$value[TAREA.class.thib.30m$class == 3])
# 
# 
# TAREA.landscape.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "landscape",
#         metric = "ta"
#     )
# 
# hist(TAREA.landscape.thib.30m$value, breaks = 30)
# 
# 
# 
# 
# ## Number of patches
# 
# 
# NUMPAT.landscape.thib.30m <-
#     sample_lsm(thibaud.2020.rast,
#                y = farm.boundary.shp,
#                level = "landscape",
#                metric = "np")
# 
# hist(NUMPAT.landscape.thib.30m$value)
# 
# 
# 
# ## Edge density
# 
# EDGE.class.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "ed"
#     )
# 
# hist(EDGE.class.thib.30m$value[EDGE.class.thib.30m$class == 3])
# 
# 
# 
# EDGE.landscape.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "landscape",
#         metric = "ed"
#     )
# 
# hist(EDGE.landscape.thib.30m$value)
# 
# 
# 
# # patch density
# 
# PD.class.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "pd"
#     )
# 
# hist(PD.class.thib.30m$value[PD.class.thib.30m$class == 3])
# 
# 
# PD.landscape.thib.30m <-
#     sample_lsm(thibaud.2020.rast,
#                y = farm.boundary.shp,
#                level = "landscape",
#                metric = "pd")
# 
# hist(PD.landscape.thib.30m$value)
# 
# 
# # percentage of landscape (PLAND)
# 
# PLAND.class.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "pland"
#     )
# 
# hist(PLAND.class.thib.30m$value[PLAND.class.thib.30m$class == 3])
# 
# 
# 
# 
# ## Diversity of land covers--Shannon's diversity
# 
# SHDI.landscape.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "landscape",
#         metric = "shdi"
#     )
# 
# hist(SHDI.landscape.thib.30m$value)
# 
# 
# 
# ## Contiguity Index
# 
# 
# CONT.class.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "contig"
#     )
# 
# hist(CONT.class.thib.30m$value[CONT.class.thib.30m$metric == "contig_mn" &
#                                  CONT.class.thib.30m$class == 3])
# 
# 
# 
# CONT.landscape.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "landscape",
#         metric = "contig"
#     )
# 
# hist(CONT.landscape.thib.30m$value[CONT.landscape.thib.30m$metric == "contig_cv"])
# 
# 
# 
# ## Aggregation Index (AI)
# 
# AI.class.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "ai"
#     )
# 
# hist(AI.class.thib.30m$value[ AI.class.thib.30m$class == 3])
# 
# 
# AI.landscape.thib.30m <-
#     sample_lsm(
#         thibaud.2020.rast,
#         y = farm.boundary.shp,
#         level = "landscape",
#         metric = "ai"
#     )
# 
# hist(AI.landscape.thib.30m$value)
# 

## Effective mesh size

# 
# EMS.class.thib.30m <-
#   sample_lsm(
#     thibaud.2020.rast,
#     y = farm.boundary.shp,
#     level = "class",
#     metric = "mesh"
#   )
# 
# hist(EMS.class.thib.30m$value[ EMS.class.thib.30m$class == 3])
# 
# 
# EMS.landscape.thib.30m <-
#   sample_lsm(
#     thibaud.2020.rast,
#     y = farm.boundary.shp,
#     level = "landscape",
#     metric = "mesh"
#   )
# 
# hist(EMS.landscape.thib.30m$value)



## ------- Looking at farm-based correlation ---------------------------

## Method using corrplot. First, build the different class tables--not all farms have all classes, and asking landscapemetrics to include all patch types causes a vector size error 

# 
# # For the forest class
# thib.forest <- tibble(
#     plot_id = FRACM.class.thib.30m$plot_id[FRACM.class.thib.30m$class == 3 &
#                                              FRACM.class.thib.30m$metric == "frac_mn"],
#     FRACM.mean.forest = FRACM.class.thib.30m$value[FRACM.class.thib.30m$class == 3 &
#                                                      FRACM.class.thib.30m$metric == "frac_mn"],
#     TAREA.forest = TAREA.class.thib.30m$value[TAREA.class.thib.30m$class ==
#                                                 3],
#     EDGE.forest = EDGE.class.thib.30m$value[EDGE.class.thib.30m$class == 3],
#     PD.forest = PD.class.thib.30m$value[PD.class.thib.30m$class == 3],
#     PLAND.forest = PLAND.class.thib.30m$value[PLAND.class.thib.30m$class == 3],
#     CONT.forest = CONT.class.thib.30m$value[CONT.class.thib.30m$metric == "contig_mn" &
#                                               CONT.class.thib.30m$class == 3],
#     AI.forest = AI.class.thib.30m$value[ AI.class.thib.30m$class == 3]
#     
# )
# 
# cor(thib.forest)
# corrplot::corrplot(cor(thib.forest[-1]), method = "shade")
# 
# 
# thib.landscape <- tibble(
#     plot_id = FRACM.landscape.thib.30m$plot_id[FRACM.landscape.thib.30m$metric=="frac_mn"],
#     
#     FRACM.mean.farm = FRACM.landscape.thib.30m$value[FRACM.landscape.thib.30m$metric=="frac_mn"],
#     TAREA.farm = TAREA.landscape.thib.30m$value,
#     NUMPAT.farm = NUMPAT.landscape.thib.30m$value,
#     PD.farm = PD.landscape.thib.30m$value,
#     EDGE.farm = EDGE.landscape.thib.30m$value,
#     SHDI.farm = SHDI.landscape.thib.30m$value,
#     CONT.mean.farm = CONT.landscape.thib.30m$value[CONT.landscape.thib.30m$metric == "contig_mn"],
#     AI.farm = AI.landscape.thib.30m$value
# )
# 
# corrplot(cor(thib.landscape[-1]), method = "shade")
# 
# 
# # Combine
# 
# all.thib <- left_join(thib.landscape, thib.forest)
# 
# all.thib[is.na(all.thib)] <- 0
# corrplot(cor(all.thib[-1]), method = "shade")
# pairs(all.thib[-1], lower.panel = panel.smooth, upper.panel = panel.cor)


## ------- Cleanup -------------------------------
# 
# if (cleanup == TRUE){
#     
#     remove(FRACM.class.thib.30m, FRACM.landscape.thib.30m,
#            AREAM.class.thib.30m, 
#            TAREA.class.thib.30m, TAREA.landscape.thib.30m,
#            NUMPAT.landscape.thib.30m,
#            EDGE.class.thib.30m, EDGE.landscape.thib.30m,
#            PD.class.thib.30m, PD.landscape.thib.30m,
#            PLAND.class.thib.30m,PLAND.landscape.thib.30m,
#            SHDI.landscape.thib.30m,
#            CONT.class.thib.30m, CONT.landscape.thib.30m, 
#            AI.class.thib.30m, AI.landscape.thib.30m, 
#            EMS.class.thib.30m, EMS.landscape.thib.30m,
#            thib.forest, thib.landscape)
#     
# }
# 



## ------- Calculate moving window metrics for Thibaud ---------------------------

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
# window.landscape.thib <- window_lsm(thibaud.2020.rast, 
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
    
    FRAC.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "mean.frac.dim.index"
        )
    
    writeRaster(FRAC.forest.window.thib, "../Xingue basin/Thibaud/FRAC.forest.5window.thib.tif")
    
    TAREA.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "total.area"
        )
    
    writeRaster(TAREA.forest.window.thib, "../Xingue basin/Thibaud/TAREA.forest.5window.thib.tif")
    
    
    NUMPAT.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "n.patches"
        )
    
    writeRaster(NUMPAT.forest.window.thib, "../Xingue basin/Thibaud/NUMPAT.forest.5window.thib.tif")
    
    EDGE.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "edge.density"
        )
    
    writeRaster(EDGE.forest.window.thib, "../Xingue basin/Thibaud/EDGE.forest.5window.thib.tif")
    
    
    PD.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "patch.density"
        )
    
    writeRaster(PD.forest.window.thib, "../Xingue basin/Thibaud/PD.forest.5window.thib.tif")
    
    PLAND.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "prop.landscape"
        )
    
    writeRaster(PLAND.forest.window.thib, "../Xingue basin/Thibaud/PLAND.forest.5window.thib.tif")
    
    COHESION.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "patch.cohesion.index"
        )
    
    writeRaster(COHESION.forest.window.thib, "../Xingue basin/Thibaud/COHESION.forest.5window.thib.tif")
    
    
    AI.forest.window.thib <-
        focal.lmetrics(
            thibaud.2020.rast,
            w = 5,
            land.value = 3,
            metric = "aggregation.index"
        )
    
    writeRaster(AI.forest.window.thib, "../Xingue basin/Thibaud/AI.forest.5window.thib.tif")
    
    
}