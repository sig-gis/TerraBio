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

## Note, for landscapemetrics package, this must be as a raster with units in meters.

all.cocoa.and.farm <- left_join(all.cocoa, all.farm) # will automatically join by plot_ID

## ------- Run PCA -------------------------------------------

## PCA assumes linear relationships between data. 

## PCA is sensitive to scaling. 

pca.all.cocoa.and.farm <- prcomp(all.cocoa.and.farm[-1], center = T, scale. = TRUE)

summary(pca.all.cocoa.and.farm)

## First Axis explains 41% of variation, second axis 13%, third 12%, fourth 8%.
## Overall, the first five have eigenvalues (sd) greater than 1, and cumulative
## proportion is 82%. With just the first two 

## Let's plot the PCA. 

plot(pca.all.cocoa.and.farm$x[,1:2], xlab="PC1 (41%)", ylab = "PC2 (13.4%)")
plot_prcomp(all.cocoa.and.farm[-1], prcomp_args = list(scale. = T, center = T))

## Use only farm metrics--since not all farms have cocoa. 

pca.all.farm <- prcomp(all.farm[-1], center = T, scale. = TRUE)
summary(pca.all.farm)

plot_prcomp(all.farm[-1], prcomp_args = list(scale. = T, center = T))

plot_prcomp(all.cocoa[-1])

kmeans.all.farm <- kmeans(x = pca.all.farm$x[,1:3], centers = 3)


plot(pca.all.farm$x[,1], pca.all.farm$x[,2], xlab="PC1 (43.6%)", ylab = "PC2 (21.5%)", col = kmeans.all.farm$cluster) # one extreme outlier for PC2, almost certainly the very large field.
plot(pca.all.farm$x[,2], pca.all.farm$x[,3], xlab="PC1 (15%)", ylab = "PC2 (11%)", col = kmeans.all.farm$cluster)

write.csv(x = pca.all.farm$rotation, file =  "pca-all-farm-rotation.csv")

farm.boundary.shp <- merge(farm.boundary.shp, tibble(farm_id = 1:150, EDGE.farm = EDGE.landscape$value, TAREA.farm = TOTALAREA$value, CONT.mn.farm = CONT.landscape$value[CONT.landscape$metric == "contig_mn"], PCA1 = pca.all.farm$x[,1], PCA2 = pca.all.farm$x[,2], PCA3 = pca.all.farm$x[,3]), by.x = "farm_id", by.y = "farm_id")

head(farm.boundary.shp)

shapefile(x = farm.boundary.shp, filename = "../Xingue basin/mu_sfx_150properties/mu_sfx_PCA_150properties.shp", overwrite = TRUE)

