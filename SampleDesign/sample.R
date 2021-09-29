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
library(RStoolbox)
library(factoextra)
library(cluster)


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
# kmeans.all.farm <- kmeans(x = all.farm[-1], centers = 3)
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

# pca.all.mb <- prcomp(all.mb[-1], center = T, scale. = TRUE)
# summary(pca.all.mb)
# 
# kmeans.all.mb <- kmeans(x = all.mb[-1], centers = 3)
# 
#     plot_prcomp(all.mb[-1], prcomp_args = list(scale. = T, center = T))
#     
#     plot(
#         pca.all.mb$x[, 1],
#         pca.all.mb$x[, 2],
#         xlab = "PC1",
#         ylab = "PC2",
#         col = kmeans.all.mb$cluster
#     ) # one extreme outlier for PC2, almost certainly the very large field.
#     plot(
#         pca.all.mb$x[, 2],
#         pca.all.mb$x[, 3],
#         xlab = "PC2",
#         ylab = "PC3",
#         col = kmeans.all.mb$cluster
#     )

## Forest metrics only
    
    pca.forest.mb <-
        prcomp(all.mb[, grepl(x = colnames(all.mb), pattern = ".forest")],
               center = T,
               scale. = TRUE)
    summary(pca.forest.mb)
    
    pca_farmsforest_mb <- tibble(
        farm_id = farm.boundary.shp$farm_id,
        PC1 = pca.forest.mb$x[ , 1],
        PC2 = pca.forest.mb$x[ , 2],
        PC3 = pca.forest.mb$x[ , 3],
        PC4 = pca.forest.mb$x[ , 4]
        )
    
    farm.boundary.shp$PC1 <- pca.forest.mb$x[ , 1]
    farm.boundary.shp$PC2 <- pca.forest.mb$x[ , 2]
    farm.boundary.shp$PC3 <- pca.forest.mb$x[ , 3]
    farm.boundary.shp$PC4 <- pca.forest.mb$x[ , 4]
    
    # write.csv(pca_farmsforest_mb,
    #           "../Xingue basin/MapBiomas/PCA/pca_farms-forest_mb.csv")
    # shapefile(
    #     farm.boundary.shp,
    #     "../Xingue basin/MapBiomas/PCA/pca_farmsforest_mb.shp",
    #     overwrite = T
    # )
    
    kmeans.forest.mb <- kmeans(x = all.mb[, grepl(x = colnames(all.mb), pattern = ".forest")], centers = 3)
    
    plot_prcomp(all.mb[, grepl(x = colnames(all.mb), pattern = ".forest")], prcomp_args = list(scale. = T, center = T))
    

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
    
## ------- Run Cluster on MapBiomaas ------------------------
    
    ## hcluster with forest
    
    m <- c( "average", "single", "complete", "ward")
    names(m) <- c( "average", "single", "complete", "ward")
    
    #function to compute agglomerative coefficient
    ac <- function(x) {
        agnes(all.mb[, grepl(x = colnames(all.mb), pattern = ".forest")], method = x, stand = TRUE)$ac
    }
    
    #calculate agglomerative coefficient for each clustering linkage method
    sapply(m, ac)
    
    ## Ward is best, but they're all very close. stand vs not using stand doesn't change. it does create chaining issues
    
    
    hcluster.forest.mb <-
        agnes(all.mb[, grepl(x = colnames(all.mb), pattern = ".forest")],
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

    
    
## ------- Run PCA on Mapbiomas MW -----------------------------------------
    
    mb.metric.stack <- stack(
        FRAC.forest.window.mb,
        TAREA.forest.window.mb,
        # EDGE.forest.window.mb,
        PLAND.forest.window.mb,
        COHESION.forest.window.mb,
        AI.forest.window.mb
    )
    
    pca.mw.mb <-
        rasterPCA(
            mb.metric.stack,
            nComp = 4,
            spca = T,
            maskCheck = F,
            filename = "../Xingue basin/MapBiomas/pca_mw_mp.tif",
            overwrite = T
        )
    summary(pca.mw.mb$model)

## ------- Run PCA on TerraClass -------------------------------------------

## Forest and farm boundary (landscape)
        
    pca.all.tc <- prcomp(all.tc[-1], center = T, scale. = TRUE)
    summary(pca.all.tc)
    
    kmeans.all.tc <- kmeans(x = all.tc[-1], centers = 3)
    
    
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
        kmeans(x = all.tc[, grepl(pattern = ".forest", x = colnames(all.tc))], centers = 3)
    
    
    plot_prcomp(all.tc[, grepl(pattern = ".forest", x = colnames(all.tc))], prcomp_args = list(scale. = T, center = T))
    
    
## -------- Run PCA on TCH ------------------------
    
## Forest and farm boundary (landscape)
    
    pca.all.tch <- prcomp(all.tch[-1], center = T, scale. = TRUE)
    summary(pca.all.tch)
    
    kmeans.all.tch <- kmeans(x = all.tch[-1], centers = 3)
    
    
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
        kmeans(x = all.tch[, grepl(pattern = ".forest", x = colnames(all.tch))], centers = 3)
    
    
    plot_prcomp(all.tch[, grepl(pattern = ".forest", x = colnames(all.tch))], prcomp_args = list(scale. = T, center = T))
    
    
    

## ------- Export ----------------------

# write.csv(x = pca.all.farm$rotation, file =  "pca-all-farm-rotation.csv")

# farm.boundary.shp <- merge(farm.boundary.shp, tibble(farm_id = 1:150, EDGE.farm = EDGE.landscape.farm.30m$value, TAREA.farm = TAREA.landscape.farm.30m$value, CONT.mn.farm = CONT.landscape.farm.30m$value[CONT.landscape$metric == "contig_mn"]), by.x = "farm_id", by.y = "farm_id")
# 
# head(farm.boundary.shp)
# head(farm.landuse.shp)

farm.boundary.shp <- sp::over(farm.landuse.shp, farm.boundary.shp, fn = NULL)

shapefile(x = farm.boundary.shp, filename = "../Xingue basin/mu_sfx_150properties/mu_sfx_PCAclust.shp", overwrite = TRUE)

