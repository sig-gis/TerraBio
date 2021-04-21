## Calculate PCA and create sampling design for TerraBio. 

##
## Import --> Metrics --> [[Sampling]]
##

## ------- Load Packages -------------------------------------

library(raster)
library(landscapemetrics)
library("dplyr")
library("maptools")
library(DataExplorer)

## ------- Load Data -----------------------------------------

#source("metrics.R")
#source("metrics_MapBiomas.R")

## Note, for landscapemetrics package, this must be as a raster with units in meters.


## ------- Run PCA on IMAFLORA -------------------------------------------

## PCA assumes linear relationships between data. 

## PCA is sensitive to scaling. 

pca.all.farm <- prcomp(all.farm[-1], center = T, scale. = TRUE)
summary(pca.all.farm)

pca.farm.forest.landscape <- prcomp(all.farm[ , 2:14], center = T, scale. = T)
summary(pca.farm.forest.landscape)

kmeans.all.farm <- kmeans(x = pca.all.farm$x[, 1:3], centers = 3)

## Let's plot the PCA. 

plot_prcomp(all.farm[-1], prcomp_args = list(scale. = T, center = T))
plot_prcomp(all.farm[, 2:14], prcomp_args = list(scale. = T, center = T))

plot(pca.all.farm$x[, 1:2], xlab = "PC1", ylab = "PC2")
plot(pca.all.farm$x[, c], xlab = "PC1", ylab = "PC2")

plot(
    pca.all.farm$x[, 1],
    pca.all.farm$x[, 2],
    xlab = "PC1",
    ylab = "PC2",
    col = kmeans.all.farm$cluster
) # one extreme outlier for PC2, almost certainly the very large field.
plot(
    pca.all.farm$x[, 2],
    pca.all.farm$x[, 3],
    xlab = "PC2",
    ylab = "PC3",
    col = kmeans.all.farm$cluster
)


## ------- Run PCA on MapBiomas -------------------------------------------

## PCA assumes linear relationships between data. 

## PCA is sensitive to scaling. 

pca.all.mb <- prcomp(all.mb[-1], center = T, scale. = TRUE)
summary(pca.all.mb)

kmeans.all.mb <- kmeans(x = pca.all.mb$x[,1:3], centers = 3)

## Let's plot the PCA. 

plot_prcomp(all.mb[-1], prcomp_args = list(scale. = T, center = T))

plot(pca.all.mb$x[, 1:2], xlab = "PC1", ylab = "PC2")
plot(pca.all.mb$x[, c(2, 3)], xlab = "PC2", ylab = "PC3")

plot(
    pca.all.mb$x[, 1],
    pca.all.mb$x[, 2],
    xlab = "PC1",
    ylab = "PC2",
    col = kmeans.all.mb$cluster
) # one extreme outlier for PC2, almost certainly the very large field.
plot(
    pca.all.mb$x[, 2],
    pca.all.mb$x[, 3],
    xlab = "PC2",
    ylab = "PC3",
    col = kmeans.all.mb$cluster
)


## ------- Run PCA on TerraClass -------------------------------------------

## PCA assumes linear relationships between data. 

## PCA is sensitive to scaling. 

pca.all.tc <- prcomp(all.tc[-1], center = T, scale. = TRUE)
summary(pca.all.tc)

kmeans.all.tc <- kmeans(x = pca.all.tc$x[,1:3], centers = 3)

## Let's plot the PCA. 

plot_prcomp(all.tc[-1], prcomp_args = list(scale. = T, center = T))


plot(pca.all.tc$x[, 1:2], xlab = "PC1", ylab = "PC2")
plot(pca.all.tc$x[, c(2, 3)], xlab = "PC2", ylab = "PC3")


plot(
    pca.all.tc$x[, 1],
    pca.all.tc$x[, 2],
    xlab = "PC1",
    ylab = "PC2",
    col = kmeans.all.tc$cluster
) # one extreme outlier for PC2, almost certainly the very large field.
plot(
    pca.all.tc$x[, 2],
    pca.all.tc$x[, 3],
    xlab = "PC2",
    ylab = "PC3",
    col = kmeans.all.tc$cluster
)




## Export

# write.csv(x = pca.all.farm$rotation, file =  "pca-all-farm-rotation.csv")

farm.boundary.shp <- merge(farm.boundary.shp, tibble(farm_id = 1:150, EDGE.farm = EDGE.landscape.farm.30m$value, TAREA.farm = TAREA.landscape.farm.30m$value, CONT.mn.farm = CONT.landscape.farm.30m$value[CONT.landscape$metric == "contig_mn"], PCA1 = pca.all.farm$x[,1], PCA2 = pca.all.farm$x[,2], PCA3 = pca.all.farm$x[,3]), by.x = "farm_id", by.y = "farm_id")

head(farm.boundary.shp)

shapefile(x = farm.boundary.shp, filename = "../Xingue basin/mu_sfx_150properties/mu_sfx_PCA_150properties.shp", overwrite = TRUE)

