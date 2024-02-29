# This analysis is based on soil samples collected in 2022 and 2023 at Cafe Apui
# and processed by EcoMol.

# This file take the eDNA species table provided by EcoMol after data processing
# and determines the TerraBio biodiversity indicators. There are also some
# species accumulation curves to verify sampling completeness.

## ----- Data ingestion & setup -----------------------------------

# libraries
library(zCompositions)
library(compositions)

# import processed data
source("CafeApui/YR2_2023/cafeapui_2023_data_processing.R")

# import functions
source("allianceBranding.R")
source("functions.R")
source("multiyear_functions.R")


## ----- Accumulation curves ----------------------------------

# Create species accumulation curves for 2022 and 2023 by land use.

# create species accumulation curves for 2022--note using silica data only.
plot(specaccum(apuiMatrixFarm2022[ grepl(pattern = "Pasture",
                                       x = rownames(apuiMatrixFarm2022)), ]),
     xlab = "Number of sites (farms)",
     ylab = "Number of ASVs",
     ylim = c(0,2000), 
     xlim = c(1,4))
plot(specaccum(apuiMatrixFarm2022[ grepl(pattern = "Forest",
                                     x = rownames(apuiMatrixFarm2022)), ]),
     add = TRUE, col = "green")
plot(specaccum(apuiMatrixFarm2022[ grepl(pattern = "SAF",
                                     x = rownames(apuiMatrixFarm2022)), ]),
     add = TRUE, col = "blue")
plot(specaccum(apuiMatrixFarm2022[ grepl(pattern = "Succession",
                                     x = rownames(apuiMatrixFarm2022)), ]),
     add = TRUE, col = "orange")

legend(x = 1, y = 1600,
       legend = c("Pasture", "Forest", "SAF", "Succession"),
       fill = c("black", "green", "blue", "orange"),
       cex = 1)


# create species accumulation curves for 2023.
plot(specaccum(apuiMatrixFarm2023[ grepl(pattern = "Pasture",
                                     x = rownames(apuiMatrixFarm2023)), ]),
     xlab = "Number of sites (farms)",
     ylab = "Number of ASVs",
     ylim = c(0,3800))
plot(specaccum(apuiMatrixFarm2023[ grepl(pattern = "Forest",
                                     x = rownames(apuiMatrixFarm2023)), ]),
     add = TRUE, col = "green")
plot(specaccum(apuiMatrixFarm2023[ grepl(pattern = "SAF",
                                     x = rownames(apuiMatrixFarm2023)), ]),
     add = TRUE, col = "blue")
plot(specaccum(apuiMatrixFarm2023[ grepl(pattern = "Succession",
                                     x = rownames(apuiMatrixFarm2023)), ]),
     add = TRUE, col = "orange")

legend(x = 1, y = 2700,
       legend = c("Pasture", "Forest", "SAF", "Succession"),
       fill = c("black", "green", "blue", "orange"),
       cex = 1)


# create species accumulation curves for 2023--one farm only
plot(specaccum(apuiMatrix2023[ grepl(pattern = "Farm 4-Pasture",
                                         x = rownames(apuiMatrix2023)), ]),
     xlab = "Number of sites (farms)",
     ylab = "Number of ASVs",
     ylim = c(0,1000))
plot(specaccum(apuiMatrix2023[ grepl(pattern = "Farm 4-Forest",
                                         x = rownames(apuiMatrix2023)), ]),
     add = TRUE, col = "green")
plot(specaccum(apuiMatrix2023[ grepl(pattern = "Farm 4-SAF",
                                         x = rownames(apuiMatrix2023)), ]),
     add = TRUE, col = "blue")
plot(specaccum(apuiMatrix2023[ grepl(pattern = "Farm 4-Succession",
                                         x = rownames(apuiMatrix2023)), ]),
     add = TRUE, col = "orange")

legend(x = 1, y = 1000,
       legend = c("Pasture", "Forest", "SAF", "Succession"),
       fill = c("black", "green", "blue", "orange"),
       cex = 1)


## ----- Proposed Indicator 1: Alpha -------------------------------------------
# PI1: Alpha diversity

#group strings to pass to different functions. 
apuiGroups2022 <- str_split_fixed(rownames(apuiMatrix2022), "-", 3)[,2]

apuiGroups2023 <- str_split_fixed(rownames(apuiMatrix2023), "-", 3)[,2]

temp1 <- str_split(rownames(apuiMatrix2022), pattern = "-", simplify = T) 
apuiFarmGroups2022 <- paste0(temp1[,1], "-", temp1[,2])
temp2 <- str_split(rownames(apuiMatrix2023), pattern = "-", simplify = T) 
apuiFarmGroups2023 <- paste0(temp2[,1], "-", temp2[,2])

remove(temp1, temp2)

# Create a table with the alpha diversity measures for each replicate. Note that
# in 2022 this analysis included rare species.

# by replicate
apuiAlphaRepl2022 <- alphaMetrics(apuiMatrix2022, 
                               groupNames = apuiGroups2022, 
                               replNames = str_sub(rownames(apuiMatrix2022), -1,-1))
apuiAlphaRepl2023 <- alphaMetrics(apuiMatrix2023,
                               groupNames = apuiGroups2023,
                               replNames = str_sub(rownames(apuiMatrix2023), -1,-1))

#grouped into farm by landuse based on replicate data
apuiAlphaFarm2022 <- alphaGroupMetrics(apuiMatrix2022,
                                       groupNames = apuiFarmGroups2022)
apuiAlphaFarm2023 <- alphaGroupMetrics(apuiMatrix2023,
                                  groupNames = apuiFarmGroups2023)

# grouped into 4 landuses based on pooled data
apuiAlphaGroup2022 <- alphaGroupMetrics(apuiMatrixFarm2022,
                                         groupNames = str_split_fixed(rownames(apuiMatrixFarm2022), "-", 2)[,2])
apuiAlphaGroup2023 <- alphaGroupMetrics(apuiMatrixFarm2023,
                                         groupNames = str_split_fixed(rownames(apuiMatrixFarm2023), "-", 2)[,2])


# 2023 in particular is very impacted by the super abundant species. 
apuiAlphaGroup2022sqrt <- alphaGroupMetrics(sqrt(apuiMatrixFarm2022),
                                             groupNames = str_split_fixed(rownames(apuiMatrixFarm2022), "-", 2)[,2])

apuiAlphaGroup2023sqrt <- alphaGroupMetrics(sqrt(apuiMatrixFarm2023),
                                             groupNames = str_split_fixed(rownames(apuiMatrixFarm2023), "-", 2)[,2])

## create graphics for Species Richness

# 2022
apuiAlphaFarm2022 %>%
    mutate(treatment = factor(str_split_fixed(siteType, "-", 2)[,2], levels = c("Pasture", "SAF", "Forest", "Succession"))) %>%
    ggplot(aes(treatment, speciesRichness)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = treatment),
        width = 0.05,
        height = 0,
        size = 4
    ) +
    labs(color = "Site Type", y = "Species Richness (2022)",
         x = "Site Type") +
    theme(legend.position = "bottom"
    ) +
    scale_y_continuous(limits = c(0, 1000)) + 
    scale_color_manual(values = supportingColorPalette) 

# 2023
apuiAlphaFarm2023 %>%
    mutate(treatment = factor(str_split_fixed(siteType, "-", 2)[,2], levels = c("Pasture", "SAF", "Forest", "Succession"))) %>%
    ggplot(aes(treatment, speciesRichness)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = treatment),
        width = 0.05,
        height = 0,
        size = 4
    ) +
    labs(color = "Site Type", y = "Species Richness (2023)",
         x = "Site Type") +
    theme(legend.position = "bottom"
    ) +
    scale_y_continuous(limits = c(0, 1000)) + 
    scale_color_manual(values = supportingColorPalette) 



## create graphics for ESR

# 2022
p1 <- apuiAlphaFarm2022 %>%
    mutate(treatment = factor(str_split_fixed(siteType, "-", 2)[,2], levels = c("Pasture", "SAF", "Forest", "Succession"))) %>%
    ggplot(aes(treatment, effectiveSR)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = treatment),
        width = 0.05,
        height = 0,
        size = 4
    ) +
    labs(color = "Site Type", y = "Effective Species Richness (2022)",
         x = "Site Type") +
    theme(legend.position = "bottom"
    ) +
    scale_y_continuous(limits = c(0, 150)) + 
    scale_color_manual(values = supportingColorPalette) 

ggsave("ApuiESR_2022.pdf",
       plot = p1,
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)

# 2023
p2 <- apuiAlphaFarm2023 %>%
    mutate(treatment = factor(str_split_fixed(siteType, "-", 2)[,2], levels = c("Pasture", "SAF", "Forest", "Succession"))) %>%
    ggplot(aes(treatment, effectiveSR)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = treatment),
        width = 0.05,
        height = 0,
        size = 4
    ) +
    labs(color = "Site Type", y = "Effective Species Richness (2023)",
         x = "Site Type") +
    theme(legend.position = "bottom"
    ) +
    scale_y_continuous(limits = c(0, 150)) + 
    scale_color_manual(values = supportingColorPalette) 

ggsave("ApuiESR_2023.pdf",
       plot = p2,
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)

grid.arrange(p1,p2,
             nrow = 1)

ggsave("ApuiESR_22_23.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)

remove(p1, p2)

# ----- Setup for Beta diversity etc. ------------------------------------------

# Create compositional matrices

#compMatrix2022 <- compMatrix(inputMatrix = apuiMatrix2022, z.warning = .99)
compFarm2022 <- compMatrix(inputMatrix = apuiMatrixFarm2022, z.warning = .99)
compType2022 <- compMatrix(inputMatrix = apuiMatrixType2022, z.warning = .99)

#compMatrix2023 <- compMatrix(inputMatrix = apuiMatrix2023, z.warning = .99)
compFarm2023 <- compMatrix(inputMatrix = apuiMatrixFarm2023, z.warning = .99)
compType2023 <- compMatrix(inputMatrix = apuiMatrixType2023, z.warning = .99)

# some columns (ASVs) removed. no columns so ok

## ----- Indicator 2: Beta diversity w/ Aitchison distance -----------------------------
# Aitchison distance uses the Euclidian distance of the compositional data that
# has been center log transformed; see
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/

## Pooled replicates for each land use and site/farm

aitchisonSite2022 <- vegdist(compFarm2022, "euc", diag = F)

min(aitchisonSite2022)
max(aitchisonSite2022)

aitchisonSite2023 <- vegdist(compFarm2023, "euc", diag = F)

min(aitchisonSite2023)
max(aitchisonSite2023)

## create plots
aitchisonSite2022 %>% aitHeatmap() 

aitchisonSite2023 %>% aitHeatmap() 


# Create the beta diversity boxplots
levelOrder = c(
    "Forest-Forest",
    "Pasture-Pasture",
    "SAF-SAF",
    "Succession-Succession",
    "Forest-Pasture",
    "Forest-SAF",
    "Forest-Succession",
    "Pasture-SAF",
    "Pasture-Succession",
    "SAF-Succession"
)

aitComparison(
    inputDist = aitchisonSite2022,
    remap = lookupSitenames2022[, c(17,17,8)], # remap needs to be old, new, type
    repeatSamples = TRUE,
    fillColor = c(
        supportingColorPalette[c(3, 1, 2, 4, 5, 6, 7, 8)],
        corporateColorPalette[c(2:3)]
    ),
    levelsPlot = levelOrder
) 


aitComparison(
    inputDist = aitchisonSite2023,
    remap = lookupSitenames2023[, c(17,17,9)], # remap needs to be old, new, type
    repeatSamples = TRUE,
    fillColor = c(
        supportingColorPalette[c(3, 1, 2, 4, 5, 6, 7, 8)],
        corporateColorPalette[c(2:3)]
    ),
    levelsPlot = levelOrder
) 



## Pooled land use types across farms
aitchisonLU2022 <- vegdist(compType2022, "euc", diag = F)

min(aitchisonLU2022)
max(aitchisonLU2022)

aitchisonLU2023 <- vegdist(compType2023, "euc", diag = F)

min(aitchisonSample2023)
max(aitchisonSample2023)

## create plots
p1 <-
    aitchisonLU2022 %>% aitHeatmap(fillColor1 = supportingColorPalette[2],
                                   fillColor2 = corporateColorPalette[4])
ggsave(
    "ApuiDistHeatmap_2022.pdf",
    plot = p1,
    device = "pdf",
    path = "OutputImages/",
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
)


p2 <- aitchisonLU2023 %>% aitHeatmap(fillColor1 = supportingColorPalette[2],
                                       fillColor2 = corporateColorPalette[4]) 
ggsave("ApuiDistHeatmap_2023.pdf",
       plot = p2,
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300
)

grid.arrange(p1, p2, nrow = 1)
ggsave("ApuiDistHeatmap_22_23.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300
)
remove(p1,p2)


## ----- Indicator 3: PCA qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# 2022 create Farm-landcover pca
pca_LCplots2022 <- compFarm2022 %>%
    PCA(., scale.unit = F, graph = F)

# 2022 PCA visualization
viz_pcaLCPlots2022 <- fviz_pca_ind(
    pca_LCplots2022,
    geom.ind = "point",
    col.ind = str_split_fixed(rownames(compFarm2022), "-", 2)[,2],
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "LC Type",
    palette = supportingColorPalette[c(3,1,2,4)]
)

ggpubr::ggpar(viz_pcaLCPlots2022,
              title = "2022 Community Composition Visualization using PCA",
              subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F
)

ggsave("ApuiPCA_2022.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)


# 2023 create Farm-landcover pca
pca_LCplots2023 <- compFarm2023 %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaLCPlots2023 <- fviz_pca_ind(
    pca_LCplots2023,
    geom.ind = "point",
    col.ind = str_split_fixed(rownames(compFarm2023), "-", 2)[,2],
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "LC Type",
    palette = supportingColorPalette[c(3,1,2,4)]
)

ggpubr::ggpar(viz_pcaLCPlots2023,
              title = "2023 Community Composition Visualization using PCA",
              subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F
)

ggsave("ApuiPCA_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)



# ----- Indicator 3: Combined 2022 and 2023 --------------------------------------------

# combine the site species matrices for 2022 and 2023. The goal is to have farm-lu rows that are renamed to include year. Then as many columns as needed, with the overlapping ASVs merged.


# fix rownames
rownames(apuiMatrixFarm2022) <-
    paste0(rownames(apuiMatrixFarm2022), "-2022")

rownames(apuiMatrixFarm2023) <-
    paste0(rownames(apuiMatrixFarm2023), "-2023")

colnames(apuiMatrixFarm2023) <- plyr::mapvalues(colnames(apuiMatrixFarm2023), from = commonASV$ASVHeader2023,
                                             to = commonASV$ASVHeader2022,
                                             warn_missing = F)

apuiMatrixJOINT <- bind_rows(apuiMatrixFarm2022, apuiMatrixFarm2023)
apuiMatrixJOINT[is.na(apuiMatrixJOINT)] <- 0


# Create a compositional matrix for the joint matrix
compMatrixJOINT <- compMatrix(inputMatrix = apuiMatrixJOINT, z.warning = .99)

# order is forest, pasture, SAF, Succession
jointcolors <- c( "#c2dd97","#8EBF3F","#fac091", "#F68B33","#fae997","#F5D226","#e3e88c", "#CAD32B")
jointcolors2 <- c( "#c2dd97","#8EBF3F","#aaaaaa", "#444444", "#fae997","#F5D226", "#bf9edf", "#8545c2")

# see here: https://www.colorhexa.com/fac091
# info on shapes for pca: https://copyprogramming.com/howto/specify-different-pointshapes-for-var-and-ind-in-fviz-pca-biplot


# create plot pcas--2022 and 2023
pca_plots <- compMatrixJOINT %>%
    PCA(., scale.unit = F, graph = F)


viz_pcaPlots <- fviz_pca_ind(
    pca_plots,
    geom.ind = "point",
    col.ind = str_split_fixed(rownames(compMatrixJOINT), "-", 2)[,2],
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    point = 16,
    mean.point = F, 
    axes = c(1,2),
    legend.title = "Site Type",
    palette = jointcolors
)
ggpubr::ggpar(viz_pcaPlots,
              title = paste0("2022 & 2023 Community Composition Visualization using PCA"),
              subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F
)



ggsave("ApuiPCA_JOINT_22_23.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)
