## Calculate metrics data for TerraBio for the sampling design using Tree Canopy
## Height data from Peter Potopov

## This script is part of a series that calculates metrics for three separate
## maps, the 150 farm data from Imaflora (date unknown), MapBiomas 2018,
## TerraClass 2014, and Tree Canopy Height 2019. The metrics are focused largely
## on forest cover as this is known to be a key landscape component for the
## butterflies and insects of interest in our eDNA sampling.

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

cleanup <- FALSE             # True will remove intermediary processing steps for cleaner environment.
write.movingwindow <- FALSE  # True will cause moving window to run & write rasters. WARNING THIS IS VERY TIME INTENSIVE


## ------- Calculate Farm-based Metrics for TCH ----------------------------

## Tree Canopy Height data can be calculated both where farms are (using
## Imaflora boundaries) and using a moving window.

## Use the sample_lsm function to sample within farm boundaries.
## https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html

## Choose height for tree cutoff

tree.cutoff = 3

## This is to recode 
isBecomes <- cbind(c(1:tree.cutoff,101, tree.cutoff:100),
                   c(rep(0,tree.cutoff),0, rep(1,101-tree.cutoff)))

tch.2019.rast <- reclassify(tch.2019.rast, rcl=isBecomes)



## Edge complexity using Fractal Dimension Index (FDI), focus on forest patches.

FRACM.class.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "frac"#, 
        #all_classes = TRUE
    )

hist(FRACM.class.tch.30m$value[FRACM.class.tch.30m$class == 1 &
                                  FRACM.class.tch.30m$metric == "frac_mn"])


FRACM.landscape.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "frac"
    )

hist(FRACM.landscape.tch.30m$value[FRACM.landscape.tch.30m$metric ==
                                       "frac_mn"])



# Total areas


TAREA.class.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "ca"
    )

hist(TAREA.class.tch.30m$value[TAREA.class.tch.30m$class == 1])


TAREA.landscape.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "ta"
    )

hist(TAREA.landscape.tch.30m$value, breaks = 30)




## Number of patches


NUMPAT.landscape.tch.30m <-
    sample_lsm(tch.2019.rast,
               y = farm.boundary.shp,
               level = "landscape",
               metric = "np")

hist(NUMPAT.landscape.tch.30m$value)



## Edge density

EDGE.class.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "ed"
    )

hist(EDGE.class.tch.30m$value[EDGE.class.tch.30m$class == 1])



EDGE.landscape.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "ed"
    )

hist(EDGE.landscape.tch.30m$value)



# patch density

PD.class.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "pd"
    )

hist(PD.class.tch.30m$value[PD.class.tch.30m$class == 1])


PD.landscape.tch.30m <-
    sample_lsm(tch.2019.rast,
               y = farm.boundary.shp,
               level = "landscape",
               metric = "pd")

hist(PD.landscape.tch.30m$value)


# percentage of landscape (PLAND)

PLAND.class.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "pland"
    )

hist(PLAND.class.tch.30m$value[PLAND.class.tch.30m$class == 1])




## Diversity of land covers--Shannon's diversity

SHDI.landscape.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "shdi"
    )

hist(SHDI.landscape.tch.30m$value)



## Contiguity Index


CONT.class.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "contig"
    )

hist(CONT.class.tch.30m$value[CONT.class.tch.30m$metric == "contig_mn" &
                                 CONT.class.tch.30m$class == 1])



CONT.landscape.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "contig"
    )

hist(CONT.landscape.tch.30m$value[CONT.landscape.tch.30m$metric == "contig_mn"])

## Aggregation Index (AI)

AI.class.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "class",
        metric = "ai"
    )

hist(AI.class.tch.30m$value[ AI.class.tch.30m$class == 1])


AI.landscape.tch.30m <-
    sample_lsm(
        tch.2019.rast,
        y = farm.boundary.shp,
        level = "landscape",
        metric = "ai"
    )

hist(AI.landscape.tch.30m$value)


## Effective mesh size
# 
# 
# EMS.class.tch.30m <-
#     sample_lsm(
#         tch.2019.rast,
#         y = farm.boundary.shp,
#         level = "class",
#         metric = "mesh"
#     )
# 
# hist(EMS.class.tch.30m$value[ EMS.class.tch.30m$class == 1])
# 
# 
# EMS.landscape.tch.30m <-
#     sample_lsm(
#         tch.2019.rast,
#         y = farm.boundary.shp,
#         level = "landscape",
#         metric = "mesh"
#     )
# 
# hist(EMS.landscape.tch.30m$value)


## ------- Looking at farm-based correlation ---------------------------

## Method using corrplot. First, build the different class tables--not all farms have all classes, and asking landscapemetrics to include all patch types causes a vector size error 


# For the forest class
tch.forest <- tibble(
    plot_id = FRACM.class.tch.30m$plot_id[FRACM.class.tch.30m$class == 1 &
                                             FRACM.class.tch.30m$metric == "frac_mn"],
    FRACM.mean.forest = FRACM.class.tch.30m$value[FRACM.class.tch.30m$class == 1 &
                                                     FRACM.class.tch.30m$metric == "frac_mn"],
    TAREA.forest = TAREA.class.tch.30m$value[TAREA.class.tch.30m$class ==
                                                1],
    EDGE.forest = EDGE.class.tch.30m$value[EDGE.class.tch.30m$class == 1],
    PD.forest = PD.class.tch.30m$value[PD.class.tch.30m$class == 1],
    PLAND.forest = PLAND.class.tch.30m$value[PLAND.class.tch.30m$class == 1],
    CONT.forest = CONT.class.tch.30m$value[CONT.class.tch.30m$metric == "contig_mn" &
                                              CONT.class.tch.30m$class == 1],
    AI.forest = AI.class.tch.30m$value[AI.class.tch.30m$class == 1]
    
)

corrplot::corrplot(cor(tch.forest[-1]), method = "shade")


tch.landscape <- tibble(
    plot_id = FRACM.landscape.tch.30m$plot_id[FRACM.landscape.tch.30m$metric=="frac_mn"],
    
    FRACM.mean.farm = FRACM.landscape.tch.30m$value[FRACM.landscape.tch.30m$metric=="frac_mn"],
    TAREA.farm = TAREA.landscape.tch.30m$value,
    NUMPAT.farm = NUMPAT.landscape.tch.30m$value,
    PD.farm = PD.landscape.tch.30m$value,
    EDGE.farm = EDGE.landscape.tch.30m$value,
    SHDI.farm = SHDI.landscape.tch.30m$value,
    CONT.mean.farm = CONT.landscape.tch.30m$value[CONT.landscape.tch.30m$metric == "contig_mn"],
    AI.farm = AI.landscape.tch.30m$value
)

corrplot(cor(tch.landscape[-1]), method = "shade")


# Combine

all.tch <- left_join(tch.landscape, tch.forest) 

all.tch[is.na(all.tch)] <- 0
corrplot(cor(all.tch[-1]), method = "shade")
pairs(all.tch[-1], lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(all.tch[ , grepl(".forest", colnames(all.tch))], lower.panel = panel.smooth, upper.panel = panel.cor)


## ------- Cleanup for farm-based -------------------------------

if (cleanup == TRUE){
    
    remove(FRACM.class.tch.30m, FRACM.landscape.tch.30m,
           AREAM.class.tch.30m, 
           TAREA.class.tch.30m, TAREA.landscape.tch.30m,
           NUMPAT.landscape.tch.30m,
           EDGE.class.tch.30m, EDGE.landscape.tch.30m,
           PD.class.tch.30m, PD.landscape.tch.30m,
           PLAND.class.tch.30m,PLAND.landscape.tch.30m,
           SHDI.landscape.tch.30m,
           CONT.class.tch.30m, CONT.landscape.tch.30m, 
           AI.class.tch.30m, AI.landscape.tch.30m, 
           EMS.class.tch.30m, EMS.landscape.tch.30m,
           tch.forest, tch.landscape)
    
}


## ------- Calculate moving window metrics for TCH ---------------------------

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
    
    FRAC.forest.window.tch <-
        focal.lmetrics(
            tch.2019.rast,
            w = 5,
            land.value = 1,
            metric = "mean.frac.dim.index"
        )
    
    writeRaster(FRAC.forest.window.tch, "../Xingue basin/TCH/FRAC.forest.5window.tch.tif")
    
     TAREA.forest.window.tch <-
        focal.lmetrics(
            tch.2019.rast,
            w = 5,
            land.value = 1,
            metric = "total.area"
        )
    
    writeRaster(TAREA.forest.window.tch, "../Xingue basin/TCH/TAREA.forest.5window.tch.tif")
    #started 8:15 PM
    
}


