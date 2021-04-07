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
pca.all.farm <- prcomp(all.farm[-1], center = T, scale. = TRUE)

summary(pca.all.cocoa.and.farm)
summary(pca.all.farm)

## First Axis explains 41% of variation, second axis 13%, third 12%, fourth 8%.
## Overall, the first five have eigenvalues (sd) greater than 1, and cumulative
## proportion is 82%. With just the first two 

## Let's plot the PCA. 

plot(pca.all.cocoa.and.farm$x[,1:2], xlab="PC1 (41%)", ylab = "PC2 (13.4%)")
plot_prcomp(all.cocoa.and.farm[-1])

plot(pca.all.farm$x[,1], pca.all.farm$x[,2], xlab="PC1 (52%)", ylab = "PC2 (15%)")
plot(pca.all.farm$x[,2], pca.all.farm$x[,3], xlab="PC1 (15%)", ylab = "PC2 (11%)")

plot_prcomp(all.farm[-1])

plot_prcomp(all.cocoa[-1])


