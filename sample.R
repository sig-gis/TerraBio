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

## PCA assumes linear relationships between data and is sensitive to scaling. 

## All metrics, including forest, tree, crop, farm boundary (landscape)
# 
# pca.all.farm <- prcomp(all.farm[-1], center = T, scale. = TRUE)
# summary(pca.all.farm)
# 
# kmeans.all.farm <- kmeans(x = pca.all.farm$x[, 1:3], centers = 3)
# 
# plot_prcomp(all.farm[-1], prcomp_args = list(scale. = T, center = T))
# 
# plot(
#     pca.all.farm$x[, 1],
#     pca.all.farm$x[, 2],
#     xlab = "PC1",
#     ylab = "PC2",
#     col = kmeans.all.farm$cluster
# ) # one extreme outlier for PC2, almost certainly the very large field.
# plot(
#     pca.all.farm$x[, 2],
#     pca.all.farm$x[, 3],
#     xlab = "PC2",
#     ylab = "PC3",
#     col = kmeans.all.farm$cluster
# )

## Only forest and farm boundary (landscape)

pca.farm.forest.landscape <- prcomp(all.farm[, grepl(pattern = ".farm|.forest", x = colnames(all.farm))], center = T, scale. = T)
summary(pca.farm.forest.landscape)

plot_prcomp(all.farm[, grepl(pattern = ".farm|.forest", x = colnames(all.farm))], prcomp_args = list(scale. = T, center = T))

## Only forest

pca.farm.forest <- prcomp(all.farm[, grepl(pattern = ".forest", x = colnames(all.farm))], center = T, scale. = T)
summary(pca.farm.forest)

plot_prcomp(all.farm[, grepl(pattern = ".forest", x = colnames(all.farm))], prcomp_args = list(scale. = T, center = T))


## ------- Run PCA on MapBiomas -------------------------------------------

## Forest and farm boundary (landscape)

pca.all.mb <- prcomp(all.mb[-1], center = T, scale. = TRUE)
summary(pca.all.mb)

kmeans.all.mb <- kmeans(x = pca.all.mb$x[,1:3], centers = 3)

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

## Forest metrics only
    
    pca.forest.mb <- prcomp(all.mb[,grepl(x = colnames(all.mb), pattern = ".forest")], center = T, scale. = TRUE)
    summary(pca.forest.mb)
    
    kmeans.forest.mb <- kmeans(x = pca.forest.mb$x[,1:3], centers = 3)
    
    plot_prcomp(all.mb[, grepl(x = colnames(all.mb), pattern = ".forest")], prcomp_args = list(scale. = T, center = T))
    
    plot(pca.forest.mb$x[, 1:2], xlab = "PC1", ylab = "PC2")
    plot(pca.forest.mb$x[, c(2, 3)], xlab = "PC2", ylab = "PC3")
    
    plot(
        pca.forest.mb$x[, 1],
        pca.forest.mb$x[, 2],
        xlab = "PC1",
        ylab = "PC2",
        col = kmeans.all.mb$cluster
    ) # one extreme outlier for PC2, almost certainly the very large field.
    plot(
        pca.forest.mb$x[, 2],
        pca.forest.mb$x[, 3],
        xlab = "PC2",
        ylab = "PC3",
        col = kmeans.all.mb$cluster
    )

## ------- Run PCA on TerraClass -------------------------------------------

## Forest and farm boundary (landscape)
        
    pca.all.tc <- prcomp(all.tc[-1], center = T, scale. = TRUE)
    summary(pca.all.tc)
    
    kmeans.all.tc <- kmeans(x = pca.all.tc$x[,1:3], centers = 3)
    
    
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

## Forest only
    
    pca.forest.tc <-
        prcomp(all.tc[, grepl(pattern = ".forest", x = colnames(all.tc))], center = T, scale. = TRUE)
    summary(pca.forest.tc)
    
    kmeans.forest.tc <-
        kmeans(x = pca.forest.tc$x[, 1:4], centers = 3)
    
    
    plot_prcomp(all.tc[, grepl(pattern = ".forest", x = colnames(all.tc))], prcomp_args = list(scale. = T, center = T))
    
    
## -------- Run PCA on TCH ------------------------
    
## Forest and farm boundary (landscape)
    
    pca.all.tch <- prcomp(all.tch[-1], center = T, scale. = TRUE)
    summary(pca.all.tch)
    
    kmeans.all.tch <- kmeans(x = pca.all.tch$x[,1:3], centers = 3)
    
    
    plot_prcomp(all.tch[-1], prcomp_args = list(scale. = T, center = T))
    
    
    plot(pca.all.tch$x[, 1:2], xlab = "PC1", ylab = "PC2")
    plot(pca.all.tch$x[, c(2, 3)], xlab = "PC2", ylab = "PC3")
    
    
    plot(
        pca.all.tch$x[, 1],
        pca.all.tch$x[, 2],
        xlab = "PC1",
        ylab = "PC2",
        col = kmeans.all.tch$cluster
    ) # one extreme outlier for PC2, almost certainly the very large field.
    plot(
        pca.all.tch$x[, 2],
        pca.all.tch$x[, 3],
        xlab = "PC2",
        ylab = "PC3",
        col = kmeans.all.tch$cluster
    )
    
    ## Forest only
    
    pca.forest.tch <-
        prcomp(all.tch[, grepl(pattern = ".forest", x = colnames(all.tch))], center = T, scale. = TRUE)
    summary(pca.forest.tch)
    
    kmeans.forest.tch <-
        kmeans(x = pca.forest.tch$x[, 1:4], centers = 3)
    
    
    plot_prcomp(all.tch[, grepl(pattern = ".forest", x = colnames(all.tch))], prcomp_args = list(scale. = T, center = T))
    
    
    

## ------- Export ----------------------

# write.csv(x = pca.all.farm$rotation, file =  "pca-all-farm-rotation.csv")

farm.boundary.shp <- merge(farm.boundary.shp, tibble(farm_id = 1:150, EDGE.farm = EDGE.landscape.farm.30m$value, TAREA.farm = TAREA.landscape.farm.30m$value, CONT.mn.farm = CONT.landscape.farm.30m$value[CONT.landscape$metric == "contig_mn"], PCA1 = pca.all.farm$x[,1], PCA2 = pca.all.farm$x[,2], PCA3 = pca.all.farm$x[,3]), by.x = "farm_id", by.y = "farm_id")

head(farm.boundary.shp)

shapefile(x = farm.boundary.shp, filename = "../Xingue basin/mu_sfx_150properties/mu_sfx_PCA_150properties.shp", overwrite = TRUE)

