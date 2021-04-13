## Calculate metrics data for TerraBio for the farm sampling design. 

## This calculates metrics for three separate maps, the 150 farm data from
## Imaflora (date unknown), MapBiomas 2018, and TerraClass 2014. The metrics are
## focused largely on forest cover as this is known to be a key landscape
## component for the butterflies and insects of interest in our eDNA sampling

##
## Import --> [[Metrics]] --> Sampling
##

## ------- Load Packages -------------------------------------

library(raster)
library(landscapemetrics)
library(corrplot)
library("dplyr")
library("maptools")

## ------- Load Data -----------------------------------------

source("import.R")

## Note, for landscapemetrics package, this must be as a raster with units in meters.


## ------- Calculate Metrics for Imaflora ----------------------------

## Imaflora data only has data where the farms are, so moving window approach
## will not work. Thus we only calculate data for each of the farms.

listlsm <- list_lsm()

## Use the sample_lsm function to sample within farm boundaries.
## https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html

## Currently this code is set up to explore the data. To operationalize further,
## see
## https://r-spatialecology.github.io/landscapemetrics/articles/getstarted.html
## and follow directions to run multiple functions at once.


# Create a tree class 

farm.tree.rast <- farm.landuse.rast[]


## Edge complexity using Fractal Dimension Index (FDI), focus on forest patches.
    
    FRACM.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "frac")
    
    hist(FRACM.class$value[FRACM.class$class==7 & FRACM.class$metric=="frac_mn"])
    

    FRACM.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "frac")

    hist(FRACM.landscape$value[FRACM.landscape$metric=="frac_mn"])


## Forest patch distribution (clustered vs. dispersed) using Aggregation index (AI)
    
    # AI.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "ai")
    # 
    # hist(AI.class$value[AI.class$class == 7])
    
    
    # AI.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "ai")
    # 
    # hist(AI.landscape$value)
    

## Mean area of forest/farms
    
    AREAM.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "area")
    
    hist(AREAM.class$value[AREAM.class$metric=="area_mn" & AREAM.class$class==7])
    
    
    # AREAM.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "area")
    # 
    # hist(AREAM.landscape$value[AREAM.landscape$metric=="area_mn"])
    


# Total areas

    CLASSAREA <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "ca")
    
    hist(CLASSAREA$value[CLASSAREA$class==7])
    
    
    TOTALAREA <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "ta")
    
    hist(TOTALAREA$value, breaks = 30)




## Number of patches

NUMPAT.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "np")

hist(NUMPAT.landscape$value)



# patch density

PD.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "pd")

hist(PD.class$value[PD.class$class == 5])



PD.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "pd")

hist(PD.landscape$value)


# percentage of landscape (PLAND)

PLAND <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "pland")

hist(PLAND$value[PLAND$class==5])



## Edge density

EDGE.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "ed")

hist(EDGE.class$value[EDGE.class$class == 5])


EDGE.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "ed")

hist(EDGE.landscape$value)



## Diversity of land covers--Shannon's diversity

SHDI <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "shdi")

hist(SHDI$value)



## Contiguity Index

CONT.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "contig")

hist(CONT.class$value[CONT.class$metric == "contig_mn" & CONT.class$class == 5])



CONT.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "contig")

hist(CONT.landscape$value[CONT.landscape$metric == "contig_cv"])






## ------- Looking at correlation ---------------------------

## Method suggested in landscapemetrics documentation; this doesn't allow for the y argument so it doesn't work.

# class.metrics <- calculate_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class")

# show_correlation(metrics, method = "pearson")


## Method using corrplot.

all.cocoa<- tibble(
    plot_id = FRACM.class$plot_id[FRACM.class$class==5 & FRACM.class$metric=="frac_mn"],
    FRACM.mean.cocoa = FRACM.class$value[FRACM.class$class==5 & FRACM.class$metric=="frac_mn"],
    AREAM.mean.cocoa = AREAM.class$value[AREAM.class$metric=="area_mn" & AREAM.class$class==5],
    CLASSAREA.cocoa = CLASSAREA$value[CLASSAREA$class==5],
    AI.cocoa = AI.class$value[AI.class$class == 5],
    
    PD.cocoa = PD.class$value[PD.class$class == 5],
    PLAND.cocoa = PLAND$value[PLAND$class==5],
    EDGE.cocoa = EDGE.class$value[EDGE.class$class == 5],
    
    CONT.cocoa = CONT.class$value[CONT.class$metric == "contig_mn" & CONT.class$class == 5]
    
)

cor(all.cocoa)
corrplot::corrplot(cor(all.cocoa[-1]), method = "shade")


all.farm <- tibble(
    plot_id = FRACM.landscape$plot_id[FRACM.landscape$metric=="frac_mn"],
    
    FRACM.mean.farm = FRACM.landscape$value[FRACM.landscape$metric=="frac_mn"],
    TOTALAREA.farm = TOTALAREA$value,
    ## AREAM.mean.farm = AREAM.landscape$value[AREAM.landscape$metric=="area_mn"],
    ## AI.farm = AI.landscape$value,
    NUMPAT.farm = NUMPAT.landscape$value,
    PD.farm = PD.landscape$value,
    EDGE.farm = EDGE.landscape$value,
    SHDI.farm = SHDI$value,
    ## ENG.farm = ENT$value,
    ## CONDENT.farm = CONDENT$value,
    CONT.mean.farm = CONT.landscape$value[CONT.landscape$metric == "contig_mn"],

)

corrplot(cor(all.farm[-1]), method = "shade")

