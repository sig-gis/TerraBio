## Import data for TerraBio prior to form analysis. Note: While waiting on
## external storage you will need to have these files downloaded with the same
## relative path location.

##
## [[Import]] --> Metrics --> Sampling
##

## ------- Load Packages -------------------------------------

library(raster, rgdal)
library(maptools)
library(dplyr) 

# panel.cor from https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/pairs

# panel.cor puts correlation in upper panels, size proportional to correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = (cex.cor) * (r+.5))
}

buffer.size <- 250

## ------- Projections ---------------------------------------

# Working in EPSG: 5880 SIRGAS 2000/Brazil Polyconic for this exercise; can be changed.

EPSG5880 <- CRS("+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Note--for the sake of time I just converted things in QGIS & added here; in
# the future can do the whole pipeline in R.




## ------- Import Data ---------------------------------------

## SHP bounding box of farms data--buffered to 1km

extent.shp <- shapefile("../Xingue basin/mu_sfx_bbox_car_boundaries/mu_sfx_bbox_5880_1kmbuf.shp")



## Detailed 150 farms data from Imaflora

farm.landuse.shp <- shapefile("../Xingue basin/mu_sfx_150properties/mu_sfx_landuse_150properties_5880.shp", verbose = T)

# farm.landuse.rast.5m <- raster("../Xingue basin/mu_sfx_150properties/mu_sfx_landuse_150properties_5880.tif")

farm.boundary.shp <- shapefile('../Xingue basin/mu_sfx_150properties/mu_sfx_boundary_150properties_5880.shp', verbose = T)
    
    farm.boundary.shp$gid <- as.integer(farm.boundary.shp$gid)
    
    farm.boundary.shp <- farm.boundary.shp[order(farm.boundary.shp$gid) , ] 
    
    # There is no number 42, which is unfortunate. Create a new ID field so that mapping is easier later.
    
    farm.boundary.shp$farm_id <- 1:150
    
    
## Create buffers for 150 farms
    
    farm.boundary.buffer <- buffer(farm.boundary.shp, width = buffer.size, dissolve = FALSE)


# Mapbiomas data

## I still can't download the 2019 Mapbiomas v 5.0 data, so using the 2018 3.0 data for now. For the pilot, I'm masking the area of study using the boundaries of the farms buffered by 1km.

mapbiomas.2018.rast <- trim(raster("../Xingue basin/MapBiomas/mapbiomas-xingue-2018.tif"))



# TerraClass data

## The most recent data available is from 2014. See https://www.terraclass.gov.br/

# terraclass.2014.rast <- trim( raster("../Xingue basin/TerraClass/TerraClass-2014-5880.tif") )


# UMD TCH data

## The most recent data available is from 2019?

tch.2019.rast <- trim( raster ("../Xingue basin/TCH/Forest_height_2019_Xingue_5880.tif"))


# Thibaud's data 

# thibaud.2020.rast <- trim(raster("../Xingue basin/Thibaud/brazil_para_new_5880.tif"))

