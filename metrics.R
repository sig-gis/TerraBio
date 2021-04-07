## Calculate metrics data for TerraBio for the form analysis. 

##
## Import --> [[Metrics]] --> Sampling
##

## ------- Load Packages -------------------------------------

library(raster, landscapemetrics, landscapetools, dplyr, maptools)

## ------- Load Data -----------------------------------------

source("import.R")

## Note, for landscapemetrics package, this must be as a raster with units in meters.

#show_patches(farm.landuse.rast)



## ------- Calculate FORM Metrics ----------------------------

listlsm <- list_lsm()

## Currently this code is set up to explore the data. To operationalize, see
## https://r-spatialecology.github.io/landscapemetrics/articles/getstarted.html
## and follow directions to run multiple functions at once.

## Use the sample_lsm function to sample within farm boundaries.
## https://r-spatialecology.github.io/landscapemetrics/reference/sample_lsm.html


## Edge complexity using Fractal Dimension Index (FDI)

FRACM.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "frac")

hist(FRACM.class$value[FRACM.class$class==5 & FRACM.class$metric=="frac_mn"])


FRACM.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "frac")

hist(FRACM.landscape$value[FRACM.landscape$metric=="frac_mn"])


## Perimeter-area ratio

PAFRAC.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "pafrac")

hist(PAFRAC.class)


## area metrics of cocoa plots/farms

AREAM.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "area")

hist(AREAM.class$value[AREAM.class$metric=="area_mn" & AREAM.class$class==5])


AREAM.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "area")

hist(AREAM.landscape$value[AREAM.landscape$metric=="area_cv"])



# Total areas

CLASSAREA <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "ca")

hist(CLASSAREA$value[CLASSAREA$class==5])


TOTALAREA <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "ta")

hist(TOTALAREA$value, breaks = 30)



## Forest patch distribution (clustered vs. dispersed) using Aggregation index (AI)

AI.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "ai")

hist(AI$value[AI$class == 5])

AI.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "ai")

hist(AI$value)


## Cohesion

COHESION.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "cohesion")

hist(COHESION.class$value[COHESION.class$class==5])

COHESION.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "cohesion")



## Number of patches

NUMPAT.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "np")

hist(NUMPAT.landscape$value)


## ------- Calculate DENSITY Metrics -------------------------

# patch density

PD.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "pd")

hist(PD.class$value[PD.class$class == 5])



PD.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "pd")

hist(PD.landscape$value)


# percentage of landscape (PLAND)

PLAND <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "pland")

hist(PLAND$value[PLAND$class==5])

## Edge density

ed


## ------- Calculate HETEROGENEITY Metrics -------------------

## Diversity of land covers--Shannon's diversity

SHDI <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "shdi")

hist(SHDI$value)

# Shannon Entropy

ENT <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "ent")

hist(ENT$value)

# Conditional Entropy

CONDENT <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "condent")

hist(CONDENT$value)

## ------- Calculate CONNECTIVITY Metrics --------------------

## Contiguity Index

CONT.class <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = "contig")

hist(CONT.class$value[CONT.class$metric == "contig_mn" & CONT.class$class == 5])



CONT.landscape <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "landscape", metric = "contig")

hist(CONT.landscape$value[CONT.landscape$metric == "contig_mn"])






## ------- Looking at correlation ---------------------------

class.metrics <- sample_lsm(farm.landuse.rast, y = farm.boundary.shp, level = "class", metric = c("contig", "pland"))

show_correlation(metrics, method = "pearson")