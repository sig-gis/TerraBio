## Calculate metrics data for TerraBio for the sampling design using Imaflora
## farm data.

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


# Run this if you don't have the landmetrics package installed.
#library(remotes)
#
#remotes::install_github("jeffreyevans/landmetrics")
library(landmetrics)

## ------- Load Data -----------------------------------------

source("import.R")

## Note, for landscapemetrics package, this must be as a raster with units in meters.

listlsm <- list_lsm()

## ------- Calculate Metrics for Imaflora ----------------------------

## Imaflora data only has data where the farms are, so moving window approach
## will not work. Thus we only calculate data for each of the farms.

## Use the sample_lsm function to sample within farm boundaries.
## https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html

## Currently this code is set up to explore the data. To operationalize further,
## see
## https://r-spatialecology.github.io/landscapemetrics/articles/getstarted.html
## and follow directions to run multiple functions at once.

# Create a tree class (12) that includes agroforestry, forest, and silvopasture. 12 = trees, 13 = ag uses, 14 = not trees/ag.

isBecomes <- cbind(c(0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11),
                   c(0, 13, 12, 14, 14, 13, 14, 12, 14, 14, 13, 12))

farm.tree.rast.5m  <- reclassify(farm.landuse.rast.5m, rcl=isBecomes)
farm.tree.rast.30m <- reclassify(farm.landuse.rast.30m, rcl=isBecomes)




## Edge complexity using Fractal Dimension Index (FDI), focus on forest patches.
    
    FRACM.class.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "frac"#, 
            #all_classes = TRUE
        )
    
    hist(FRACM.class.farm.30m$value[FRACM.class.farm.30m$class == 7 &
                                        FRACM.class.farm.30m$metric == "frac_mn"])
    
    FRACM.class.tree.30m <-
        sample_lsm(
            farm.tree.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "frac"
        )
    
    hist(FRACM.class.tree.30m$value[FRACM.class.tree.30m$class == 12 &
                                        FRACM.class.tree.30m$metric == "frac_mn"])
    
    
    
    FRACM.landscape.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "landscape",
            metric = "frac"
        )
    
    hist(FRACM.landscape.farm.30m$value[FRACM.landscape.farm.30m$metric ==
                                            "frac_mn"])



# Total areas
    
    
    TAREA.class.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "ca"
        )
    
    hist(TAREA.class.farm.30m$value[TAREA.class.farm.30m$class == 7])
    
    TAREA.class.tree.30m <-
        sample_lsm(
            farm.tree.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "ca"
        )
    
    hist(TAREA.class.tree.30m$value[TAREA.class.tree.30m$class == 12])
    hist(TAREA.class.tree.30m$value[TAREA.class.tree.30m$class == 13])
    
    
    TAREA.landscape.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "landscape",
            metric = "ta"
        )
    
    hist(TAREA.landscape.farm.30m$value, breaks = 30)
    



## Number of patches
    
    
    NUMPAT.landscape.farm.30m <-
        sample_lsm(farm.landuse.rast.30m,
                   y = farm.boundary.shp,
                   level = "landscape",
                   metric = "np")
    
    hist(NUMPAT.landscape.farm.30m$value)


    
## Edge density

    EDGE.class.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "ed"
        )
    
    hist(EDGE.class.farm.30m$value[EDGE.class.farm.30m$class == 7])
    
    EDGE.class.tree.30m <-
        sample_lsm(
            farm.tree.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "ed"
        )
    
    hist(EDGE.class.tree.30m$value[EDGE.class.tree.30m$class == 12])
    
    EDGE.landscape.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "landscape",
            metric = "ed"
        )
    
    hist(EDGE.landscape.farm.30m$value)



# patch density

    PD.class.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "pd"
        )
    
    hist(PD.class.farm.30m$value[PD.class.farm.30m$class == 7])
    
    PD.class.tree.30m <-
        sample_lsm(
            farm.tree.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "pd"
        )
    
    hist(PD.class.tree.30m$value[PD.class.tree.30m$class == 12])
    
    
    
    PD.landscape.farm.30m <-
        sample_lsm(farm.landuse.rast.30m,
                   y = farm.boundary.shp,
                   level = "landscape",
                   metric = "pd")
    
    hist(PD.landscape.farm.30m$value)


# percentage of landscape (PLAND)

    PLAND.class.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "pland"
        )
    
    hist(PLAND.class.farm.30m$value[PLAND.class.farm.30m$class == 7])
    
    PLAND.class.tree.30m <-
        sample_lsm(
            farm.tree.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "pland"
        )
    
    hist(PLAND.class.tree.30m$value[PLAND.class.tree.30m$class == 12])
    

## Diversity of land covers--Shannon's diversity

    SHDI.landscape.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "landscape",
            metric = "shdi"
        )
    
    hist(SHDI.landscape.farm.30m$value)



## Contiguity Index
    
    
    CONT.class.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "contig"
        )
    
    hist(CONT.class.farm.30m$value[CONT.class.farm.30m$metric == "contig_mn" &
                                       CONT.class.farm.30m$class == 7])
    
    
    CONT.class.tree.30m <-
        sample_lsm(
            farm.tree.rast.30m,
            y = farm.boundary.shp,
            level = "class",
            metric = "contig"
        )
    
    hist(CONT.class.tree.30m$value[CONT.class.tree.30m$metric == "contig_mn" &
                                       CONT.class.tree.30m$class == 12])
    hist(CONT.class.tree.30m$value[CONT.class.tree.30m$metric == "contig_mn" &
                                       CONT.class.tree.30m$class == 13])    
    
    
    
    CONT.landscape.farm.30m <-
        sample_lsm(
            farm.landuse.rast.30m,
            y = farm.boundary.shp,
            level = "landscape",
            metric = "contig"
        )
    
    hist(CONT.landscape.farm.30m$value[CONT.landscape.farm.30m$metric == "contig_cv"])
    

    

## ------- Looking at correlation for Imaflora ---------------------------

## Method suggested in landscapemetrics documentation; this doesn't allow for the y argument so it doesn't work.

# class.metrics <- calculate_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class")

# show_correlation(metrics, method = "pearson")


## Method using corrplot. First, build the different class tables--not all farms have all classes, and asking landscapemetrics to include all patch types causes a vector size error (36--76GB, so can run if available memory is larger.)

    
# For the forest class
    farm.forest <- tibble(
        plot_id = FRACM.class.farm.30m$plot_id[FRACM.class.farm.30m$class == 7 &
                                                   FRACM.class.farm.30m$metric == "frac_mn"],
        FRACM.mean.forest = FRACM.class.farm.30m$value[FRACM.class.farm.30m$class == 7 &
                                                           FRACM.class.farm.30m$metric == "frac_mn"],
        TAREA.forest = TAREA.class.farm.30m$value[TAREA.class.farm.30m$class ==
                                                      7],
        EDGE.forest = EDGE.class.farm.30m$value[EDGE.class.farm.30m$class == 7],
        PD.forest = PD.class.farm.30m$value[PD.class.farm.30m$class == 7],
        PLAND.forest = PLAND.class.farm.30m$value[PLAND.class.farm.30m$class == 7],
        CONT.forest = CONT.class.farm.30m$value[CONT.class.farm.30m$metric == "contig_mn" &
                                                    CONT.class.farm.30m$class == 7]
        
    )

    cor(farm.forest)
    corrplot::corrplot(cor(farm.forest[-1]), method = "shade")

## For the tree class

    farm.tree <- tibble(
        plot_id = FRACM.class.tree.30m$plot_id[FRACM.class.tree.30m$class == 12 &
                                                   FRACM.class.tree.30m$metric == "frac_mn"],
        FRACM.mean.tree = FRACM.class.tree.30m$value[FRACM.class.tree.30m$class == 12 &
                                                           FRACM.class.tree.30m$metric == "frac_mn"],
        TAREA.tree = TAREA.class.tree.30m$value[TAREA.class.tree.30m$class ==
                                                      12],
        EDGE.tree = EDGE.class.tree.30m$value[EDGE.class.tree.30m$class == 12],
        PD.tree = PD.class.tree.30m$value[PD.class.tree.30m$class == 12],
        PLAND.tree = PLAND.class.tree.30m$value[PLAND.class.tree.30m$class == 12],
        CONT.tree = CONT.class.tree.30m$value[CONT.class.tree.30m$metric == "contig_mn" &
                                                    CONT.class.tree.30m$class == 12]
        
    )
    
    corrplot::corrplot(cor(farm.tree[-1]), method = "shade")
    # correlation relationships for forest & trees look v. close--

## For the crop class
    
    farm.crop <- tibble(
        plot_id = FRACM.class.tree.30m$plot_id[FRACM.class.tree.30m$class == 13 &
                                                   FRACM.class.tree.30m$metric == "frac_mn"],
        FRACM.mean.crop = FRACM.class.tree.30m$value[FRACM.class.tree.30m$class == 13 &
                                                         FRACM.class.tree.30m$metric == "frac_mn"],
        TAREA.crop = TAREA.class.tree.30m$value[TAREA.class.tree.30m$class ==
                                                    13],
        EDGE.crop = EDGE.class.tree.30m$value[EDGE.class.tree.30m$class == 13],
        PD.crop = PD.class.tree.30m$value[PD.class.tree.30m$class == 13],
        PLAND.crop = PLAND.class.tree.30m$value[PLAND.class.tree.30m$class == 13],
        CONT.crop = CONT.class.tree.30m$value[CONT.class.tree.30m$metric == "contig_mn" &
                                                  CONT.class.tree.30m$class == 13]
        
    )
    
    corrplot::corrplot(cor(farm.crop[-1]), method = "shade")
    # correlation relationships for forest & trees & crops look v. close--my
    # guess is that we can stick to just forest. Forest and crop likely inverse
    # relationship, with area of agroforestry maybe not enough to matter.
    




farm.landscape <- tibble(
    plot_id = FRACM.landscape.farm.30m$plot_id[FRACM.landscape.farm.30m$metric=="frac_mn"],
    
    FRACM.mean.farm = FRACM.landscape.farm.30m$value[FRACM.landscape.farm.30m$metric=="frac_mn"],
    TAREA.farm = TAREA.landscape.farm.30m$value,
    NUMPAT.farm = NUMPAT.landscape.farm.30m$value,
    PD.farm = PD.landscape.farm.30m$value,
    EDGE.farm = EDGE.landscape.farm.30m$value,
    SHDI.farm = SHDI.landscape.farm.30m$value,
    CONT.mean.farm = CONT.landscape.farm.30m$value[CONT.landscape.farm.30m$metric == "contig_mn"],

)

corrplot(cor(farm.landscape[-1]), method = "shade")


# Combine into one table

all.farm <- left_join(farm.landscape, farm.forest) %>%
    left_join(., farm.tree) %>%
    left_join(., farm.crop)# will automatically join by plot_ID

all.farm[is.na(all.farm)] <- 0
corrplot(cor(all.farm[-1]), method = "shade")
corrplot(cor(all.farm[ , 2:14]), method = "shade")
pairs(all.farm[ , 2:14], lower.panel = panel.smooth, upper.panel = panel.cor)


## ------- Cleanup -------------------------------

remove(FRACM.class.farm.30m, FRACM.class.tree.30m, FRACM.landscape.farm.30m,
       AREAM.class.farm.30m, 
       TAREA.class.farm.30m, TAREA.class.tree.30m, TAREA.landscape.farm.30m,
       NUMPAT.landscape.farm.30m,
       EDGE.class.farm.30m, EDGE.class.tree.30m, EDGE.landscape.farm.30m,
       PD.class.farm.30m, PD.class.tree.30m, PD.landscape.farm.30m,
       PLAND.class.farm.30m, PLAND.class.tree.30m, PLAND.landscape.farm.30m,
       SHDI.landscape.farm.30m,
       CONT.class.farm.30m, CONT.class.tree.30m, CONT.landscape.farm.30m)

