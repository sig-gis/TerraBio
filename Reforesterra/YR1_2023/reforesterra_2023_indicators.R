# This analysis is based on soil samples collected at Reforesterra and processed by
# EcoMol. There are three main land uses: CF = counter factual, R =
# control/forest, I = intervention; and some variants: CF2 = riparian
# counterfatual, IA and IB + ??.  The 2023 analysis is the first year.

# This analysis includes data from 2023

# This file take the eDNA species table provided by EcoMol and determines the
# TerraBio biodiversity indicators. There are also some species accumulation
# curves to verify sampling completeness.

# Code written Dec. 2023 by Karen Dyson



## ----- Data ingestion & setup -----------------------------------

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(iNEXT)

# Ingest codes
source("allianceBranding.R")
source("functions.R")
source("multiyear_functions.R")
source("../../../RCode/R_Scripts/triplet_fixer.R") # From my github repository
source("../../../RCode/R_Scripts/repeat_multipatt.R") # ditto

source("Reforesterra/YR1_2023/reforesterra_2023_data_processing.R")

## ----- Accumulation curves ----------------------------------

# create species accumulation curves
plot(specaccum(rterra2023Matrix[ grepl(pattern = "I",
                                      x = rownames(rterra2023Matrix)), ]),
     xlab = "Number of replicates 2022",
     ylab = "Number of species",
     ylim = c(0,1000))
plot(specaccum(rterra2023Matrix[ grepl(pattern = "-R-",
                                      x = rownames(rterra2023Matrix)), ]),
     add = TRUE, col = "green")
plot(specaccum(rterra2023Matrix[ grepl(pattern = "CF-",
                                      x = rownames(rterra2023Matrix)), ]),
     add = TRUE, col = "blue")
plot(specaccum(rterra2023Matrix[ grepl(pattern = "CF2",
                                      x = rownames(rterra2023Matrix)), ]),
     add = TRUE, col = "purple")

legend(x = 1, y = 700,
       legend = c("Intervention", "Forest", "Pasture", "Riparian Pasture"),
       fill = c("black", "green", "blue", "purple"),
       cex = 1)



## ----- Proposed Indicator 1: Alpha -------------------------------------------
# PI1: Alpha diversity

lookupSitenames2023$sample.code[!(lookupSitenames2023$sample.code %in% removedSites) & !(lookupSitenames2023$sample.code %in% rownames(rterra2023Matrix))]


# Create a table with the alpha diversity measures for each replicate. Note that
# in 2022 this analysis included rare species.

rterrAlpha2023 <- alphaMetrics(
    rterra2023Matrix,
    groupNames = str_split_fixed(rownames(rterra2023Matrix), "-", 5)[, 4],
    replNames = str_split_fixed(rownames(rterra2023Matrix), "-", 5)[, 5]
) %>%
    mutate(siteType = factor(
        siteType,
        levels = c("CF", "CF2", "I", "IA", "IB", "R"),
        labels = c(
            "Counterfactual",
            "Riparian CF",
            "Intervention",
            "Interv. A",
            "Interv. B",
            "Forest"
        )
    ))

rterrAlpha2023$siteNumber <- str_split_fixed(rterrAlpha2023$siteNames, "-",5)[,3]

    


rterraAlphaGroup2023 <- alphaGroupMetrics(rterraSite2023Matrix,
                                          groupNames = str_split_fixed(rownames(rterraSite2023Matrix), "-", 2)[, 1]) %>%
    mutate(siteType = factor(
        siteType,
        levels = c("CF", "CF2", "I", "IA", "IB", "R"),
        labels = c(
            "Counterfactual",
            "Riparian CF",
            "Intervention",
            "Interv. A",
            "Interv. B",
            "Forest"
        )
    ))

rterraAlphaGroup2023sqrt <-
    alphaGroupMetrics(sqrt(rterraSite2023Matrix),
                      groupNames = str_split_fixed(rownames(rterraSite2023Matrix), "-", 2)[, 1]) %>%
    mutate(siteType = factor(
        siteType,
        levels = c("CF", "CF2", "I", "IA", "IB", "R"),
        labels = c(
            "Counterfactual",
            "Riparian CF",
            "Intervention",
            "Interv. A",
            "Interv. B",
            "Forest"
        )
    ))




# Graph species richness

rterrAlpha2023 %>%
    ggplot() +
    geom_boxplot(aes(siteType, speciesRichness),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, speciesRichness, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Site Type", y = "Raw Species Richness 2023",
         x = "Site Type") +
    theme(legend.position="bottom") +
    scale_color_manual(values = supportingColorPalette)


rterrAlpha2023 %>%
    ggplot() +
    geom_boxplot(aes(siteType, speciesRichness),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, speciesRichness, color = siteType),
        width = 0.1,
        height = 0
    ) +
    facet_grid(cols = vars(siteNumber))+
    labs(color = "Site Type", y = "Raw Species Richness 2023",
         x = "Site Type") +
    theme(legend.position="bottom", 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    scale_color_manual(values = supportingColorPalette)

# graph for group SR
rterraAlphaGroup2023 %>% 
    ggplot(
        aes(siteType, speciesRichness, fill = siteType)
    ) +
    geom_bar(stat = "identity") +
    labs(fill = "Site Type", 
         y = "Raw Species Richness 2023",
         x = "Site Type") +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = supportingColorPalette) +
    geom_text(stat='identity', aes(label=round(speciesRichness)),position = position_stack(vjust = 0.5))




# Graph for ESR

rterrAlpha2023 %>%
    ggplot() +
    geom_boxplot(aes(siteType, effectiveSR),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, effectiveSR, color = siteType),
        width = 0.05,
        height = 0, size = 4
    ) +
    labs(color = "Site Type", y = "Effective Species Richness",
         x = "Site Type") + 
    theme(legend.position="bottom")+
    scale_color_manual(values = supportingColorPalette)

# graph for group ESR
rterraAlphaGroup2023 %>% 
    ggplot(
        aes(siteType, effectiveSR, fill = siteType)
    ) +
    geom_bar(stat = "identity") +
    labs(fill = "Site Type", 
         y = "Effective Species Richness 2023",
         x = "Site Type") +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = supportingColorPalette) +
    geom_text(stat='identity', aes(label=round(effectiveSR)),position = position_stack(vjust = 0.5))

rterraAlphaGroup2023sqrt %>% 
    ggplot(
        aes(siteType, effectiveSR, fill = siteType)
    ) +
    geom_bar(stat = "identity") +
    labs(fill = "Site Type", 
         y = "Effective Species Richness 2023",
         x = "Site Type") +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = supportingColorPalette) +
    geom_text(stat='identity', aes(label=round(effectiveSR)),position = position_stack(vjust = 0.5))

ggsave("RterraESR_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)



# ----- Setup for Beta diversity etc. ------------------------------------------

# Create compositional matrices

compMatrix2023 <- compMatrix(inputMatrix = rterra2023Matrix, z.warning = .999)
    # row 11, 13, 14, 15, 19, 20, 22 removed with .99
compSite2023 <- compMatrix(inputMatrix = rterraSite2023Matrix, z.warning = .999)

compType2023 <- compMatrix(inputMatrix = rterr)
# complete this--need matrix grouped by land cover

# updated package is more aggressive about throwing away data--discourage it with the z.warning.





## ----- Indicator 2: Beta diversity w/ Aitchison distance -----------------------------
# Aitchison distance uses the Euclidian distance of the compositional data that
# has been center log transformed; see
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/

# For each replicate--internal use only, really.
aitchisonReplicate2023 <- vegdist(compMatrix2023, "euc", diag = F)

min(aitchisonReplicate2023)
max(aitchisonReplicate2023)

aitchisonReplicate2023 %>% aitHeatmap()

levelOrder = c(
    "Counterfactual-Counterfactual",
    "Forest-Forest",
    "Restoration-Restoration",
    "Syntropic-Syntropic",
    "Counterfactual-Forest",
    "Counterfactual-Restoration",
    "Counterfactual-Syntropic",
    "Forest-Restoration",
    "Forest-Syntropic",
    "Restoration-Syntropic"
)
aitComparison(
    inputDist = aitchisonReplicate2023,
    remap = lookupSitenames2023[, c(1,6,4)], # remap needs to be old, new, type
    repeatSamples = TRUE,
    fillColor = c(
        supportingColorPalette[c(1, 5, 2, 6, 7, 3, 8)],
        corporateColorPalette[c(1, 2)],
        supportingColorPalette[4]
    ),
    levelsPlot = levelOrder
) # same basic pattern sqrt vs not





## For each sample (pooled replicates by land use type)

aitchisonSample2023 <- vegdist(compType2023, "euc", diag = F)

min(aitchisonSample2023)
max(aitchisonSample2023)

## create plots

treatmentHeatmap <-
    aitchisonSample2023 %>% aitHeatmap(fillColor1 = supportingColorPalette[2],
                                      fillColor2 = corporateColorPalette[4]) 


ggsave("HortaDistHeatmap_2023.pdf",
       plot = treatmentHeatmap,
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300
)

remove(treatmentHeatmap)


## ----- Indicator 3: PCA qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/


# create plot pcas--2023 only
pca_plots <- compMatrix2023 %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaPlots <- fviz_pca_ind(
    pca_plots,
    geom.ind = "point",
    col.ind = hortaGroups2023,
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "Site Type",
    palette = supportingColorPalette[c(1,3,4,2)]
)
ggpubr::ggpar(viz_pcaPlots,
              title = paste0("Community Composition Visualization using PCA"),#, hortaSubset),
              subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F, orientation = "horizontal"
)
ggsave("HortaPCA_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)


viz_pcaPlots_contrib <- fviz_contrib(pca_plots, choice = "ind", axes = 1:2)

# fviz_pca_biplot(pca_plots,
#                 # Sites
#                 col.ind = hortaGroups2023,
#                 addEllipses = T,
#                 ellipse.type = "convex",
#                 label = "var",
#                 repel = T,
#                 max.overlaps = 5,
#                 alpha.var ="contrib", )

remove(pca_plots, viz_pcaPlots, viz_pcaPlots_contrib)

# ----- Indicator 3: Combined 2022 and 2023 --------------------------------------------

# combine the site species matrices for 2022 and 2023. The goal is to have 23 rows that are renamed to include year. Then as many columns as neeeded, with the overlapping ASVs merged.


# fix rownames--warning if you switch to silica
rownames(hortaMatrix2022) <-
    str_sub(rownames(hortaMatrix2022),-3, -1) %>%
    paste0(., "-2022")
rownames(hortaMatrix2023) <-
    str_sub(rownames(hortaMatrix2023),-2, -1) %>%
    paste0(., "-2023")
    
colnames(hortaMatrix2023) <- plyr::mapvalues(colnames(hortaMatrix2023), from = commonASV$ASVHeader2023,
                                            to = commonASV$ASVHeader2022,
                                            warn_missing = F)

hortaMatrixJOINT <- bind_rows(hortaMatrix2022, hortaMatrix2023)
hortaMatrixJOINT[is.na(hortaMatrixJOINT)] <- 0


# Create a compositional matrix for the joint matrix
compMatrixJOINT <- compMatrix(inputMatrix = hortaMatrixJOINT, z.warning = .99)

# order is counterfactual, forest, restoration, syntropic
jointcolors <- c("#fac091", "#F68B33", "#c2dd97","#8EBF3F","#e3e88c", "#CAD32B","#fae997","#F5D226")
jointcolors2 <- c( "#aaaaaa", "#444444", "#c2dd97","#8EBF3F","#fae997","#F5D226", "#bf9edf", "#8545c2")

# see here: https://www.colorhexa.com/fac091
# info on shapes for pca: https://copyprogramming.com/howto/specify-different-pointshapes-for-var-and-ind-in-fviz-pca-biplot


# create plot pcas--2022 and 2023
    pca_plots <- compMatrixJOINT %>%
        PCA(., scale.unit = F, graph = F)
    
    viz_pcaPlots <- fviz_pca_ind(
        pca_plots,
        geom.ind = "point",
        col.ind = c(paste0(hortaGroups2022,"-2022"), paste0(hortaGroups2023,"-2023")),
        addEllipses = T,
        ellipse.type = "confidence",
        pointsize = 3,
        point = 16,
        mean.point = F,
        legend.title = "Site Type",
        palette = jointcolors2
    )
    ggpubr::ggpar(viz_pcaPlots,
                  title = paste0("Community Composition Visualization using PCA"),#, hortaSubset),
                  subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F
    )
    ggsave("HortaPCA_JOINT.pdf",
           plot = last_plot(),
           device = "pdf",
           path = "OutputImages/",
           width = 8,
           height = 5,
           units = "in",
           dpi = 300
    )
    