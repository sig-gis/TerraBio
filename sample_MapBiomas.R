## Calculate PCA and create sampling design for TerraBio. 

##
## Import --> Metrics --> [[Sampling]]
##

## ------- Load Packages -------------------------------------

library(raster)
library(landscapemetrics)
library("maptools")
library(DataExplorer)
library(RStoolbox)
library(factoextra)
library(cluster)
library(geosphere)
library("dplyr")
library(sf)
library(rgeos)


## ------- Load Data ------------------------------------------------------

source("metrics_MapBiomas.R")

## Note, for landscapemetrics package, you must have a raster with units in meters.


## ------- Identify shaded cocoa ------------------------------------------

# Mean height for each field

mean.height <- extract(tch.2019.rast, farm.landuse.shp, fun = mean, na.rm = TRUE)

farm.landuse.shp$mean.ht <- mean.height[,1]

farm.landuse.shp$shd.cocoa <-
    ifelse(
        farm.landuse.shp$mean.ht > 7.4 &
            farm.landuse.shp$landuse_cl == "cocoa",
        "shaded",
        ifelse(
            farm.landuse.shp$mean.ht <= 7.4 &
                farm.landuse.shp$landuse_cl == "cocoa"
            ,
            "sun",
            "notcocoa"
        )
    )

table(farm.landuse.shp$shd.cocoa)


## and now specify which farms have shaded cocoa

farm.boundary.shp$shd.cocoa <- ifelse(farm.boundary.shp$id_prod %in%
                                          farm.landuse.shp$id_prod[farm.landuse.shp$shd.cocoa == "shaded"],
                                      yes = "shaded",
                                      no = "notshaded")

## ------- Run PCA on MapBiomas -------------------------------------------

## Forest metrics from the farm boundaries & buffers

pca.forest.mb <-
    prcomp(all.mb[-1],
           center = T,
           scale. = TRUE)
summary(pca.forest.mb)

plot_prcomp(all.mb[-1], prcomp_args = list(scale. = T, center = T))


## Now we can remove the variables that do not appear strong in the first two axes. This includes NUMPAT.buffer, FRACM, and EDGE from the farm boundaries. 

reduced.mb <- all.mb[ , !(grepl(x = colnames(all.mb), pattern = "NUMPAT.buffer|FRACM|EDGE|PD.forest")) ]

pca.forest.mb <-
    prcomp(reduced.mb[-1],
           center = T,
           scale. = TRUE)
summary(pca.forest.mb)

plot_prcomp(reduced.mb[-1], prcomp_args = list(scale. = T, center = T))

# very nice!



## ------- Run Cluster on MapBiomas ------------------------

## hcluster with forest

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
    agnes(reduced.mb[-1], method = x, stand = TRUE)$ac
}

#calculate agglomerative coefficient for each clustering linkage method

sapply(m, ac)

## Ward is best, but they're all very close. stand vs not using stand doesn't change. it does create chaining issues


hcluster.forest.mb <-
    agnes(reduced.mb[-1],
          method = "ward", stand = TRUE)
plot(hcluster.forest.mb)
groups <- cutree(hcluster.forest.mb, k=5)

farm.boundary.shp$group <- groups

plot(
    pca.forest.mb$x[, 1],
    pca.forest.mb$x[, 2],
    xlab = "PC1",
    ylab = "PC2",
    col = groups
) # one extreme outlier for PC2, almost certainly the very large field.
plot(
    pca.forest.mb$x[, 2],
    pca.forest.mb$x[, 3],
    xlab = "PC2",
    ylab = "PC3",
    col = groups
)

fviz_dend(hcluster.forest.mb, k = 5)



## ------- Run summary statistics -----------------------------------------

## Find distance to river

xingu <- shapefile("../Xingue basin/Xingu_river_5880.shp")

temp <- t(rgeos::gDistance(farm.boundary.shp, xingu, byid = TRUE))
attr(temp, "dimnames") <- NULL

farm.boundary.shp$xingu.dist <- temp[,1]

farm.boundary.df <- as(farm.boundary.shp, "data.frame")

group_by(farm.boundary.df, shd.cocoa) %>% dplyr::summarise(
    count = length(id_prod),
    count.unique = length(unique(id_prod)),
    mean.area = mean(area),
    sd.area = sd(area),
    mean.dist = mean(xingu.dist),
    sd.dist = sd(xingu.dist)
)

buffer.metric.summary <-
    left_join(farm.boundary.df, all.mb, c("farm_id" = "plot_id")) %>%
    group_by(., shd.cocoa) %>%
    select(contains("buffer")) %>%
    dplyr::summarise_all(.funs = c("mean", "sd"))

farm.metric.summary <-
    left_join(farm.boundary.df, all.mb, c("farm_id" = "plot_id")) %>%
    group_by(., shd.cocoa) %>%
    select(contains("forest")) %>%
    dplyr::summarise_all(.funs = c("mean", "sd"))

head(all.mb)


## ------- Write files ------------------------------------------

farm.boundary.shp$area_sqm <- area(farm.boundary.shp)
farm.boundary.shp$area <- area(farm.boundary.shp) / 10000

farm.landuse.shp$area_ha <- area(farm.landuse.shp) / 10000



## conducting sampling right now in QGIS b.c it is faster.


shapefile(
    farm.boundary.shp,
    "../Xingue basin/MapBiomas/R_output/mu_sfx_boundary_clust_5880.shp",
    overwrite = T
)

shapefile(
    farm.landuse.shp,
    "../Xingue basin/MapBiomas/R_output/mu_sfx_landuse_clust_5880.shp",
    overwrite = T
)
