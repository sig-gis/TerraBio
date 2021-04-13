## Import data for TerraBio prior to form analysis. Note: While waiting on
## external storage you will need to have these files downloaded with the same
## relative path location.

##
## [[Import]] --> Metrics --> Sampling
##

## ------- Load Packages -------------------------------------

library(raster, rgdal)
library(maptools)
library("dplyr") # dplyr wants to be in quotes, no idea why.

## ------- Projections ---------------------------------------

# Working in EPSG: 5880 SIRGAS 2000/Brazil Polyconic for this exercise; can be changed.

EPSG5880 <- CRS("+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Note--for the sake of time I just converted things in QGIS & added here; in
# the future can do the whole pipeline in R.




## ------- Import Data ---------------------------------------

## SHP bounding box of farms data--buffered to 1km

extent.shp <- shapefile("../Xingue basin/mu_sfx_bbox_car_boundaries/mu_sfx_bbox_5880_1kmbuf.shp")

## Rasterized bounding box of farms data--buffered to 1km at 5m resolution

extent.rast <- raster("../Xingue basin/mu_sfx_bbox_car_boundaries/mu_sfx_bbox_5880_1kmbuf.tif")



## Detailed 150 farms data from Imaflora

farm.landuse.shp <- shapefile("../Xingue basin/mu_sfx_150properties/mu_sfx_landuse_150properties_5880.shp", verbose = T)

farm.landuse.rast.5m <- raster("../Xingue basin/mu_sfx_150properties/mu_sfx_landuse_150properties_5880.tif")
farm.landuse.rast.30m <- raster("../Xingue basin/mu_sfx_150properties/mu_sfx_landuse_150properties_5880_30m.tif")
farm.tree.rast.30m

farm.boundary.shp <- shapefile('../Xingue basin/mu_sfx_150properties/mu_sfx_boundary_150properties_5880.shp', verbose = T)


farm.boundary.shp$gid <- as.integer(farm.boundary.shp$gid)

farm.boundary.shp <- farm.boundary.shp[order(farm.boundary.shp$gid) , ] 

# There is no number 42, which is unfortunate. Create a new ID field so that mapping is easier later.

farm.boundary.shp$farm_id <- 1:150


# Mapbiomas data

## I still can't download the 2019 Mapbiomas v 5.0 data, so using the 2018 3.0 data for now. For the pilot, I'm masking the area of study using the boundaries of the farms buffered by 1km.

mapbiomas.2018.rast <- raster("../Xingue basin/MapBiomas/mapbiomas-xingue-2018.tif") 



# TerraClass data

## The most recent data available is from 2014

terraclass.2014.rast <- raster("../Xingue basin/TerraClass/TerraClass-2014-5880.tif")
