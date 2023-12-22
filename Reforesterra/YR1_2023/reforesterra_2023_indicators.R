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

compType2023 <- compMatrix(inputMatrix = rterraLC2023Matrix, z.warning = .99)

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
aitchisonReplicate2023

levelOrder = c(
    "CF-CF",
    "CF-CF2",
    "CF-I",
    "CF-IA",
    "CF-R",
    "CF2-CF2",
    "CF2-I",
    "CF2-IA",
    "CF2-R",
    "I-I",
    "I-IA",
    "I-R",
    "IA-IA",
    "IA-R",
    "R-R"
)
aitComparison(
    inputDist = aitchisonReplicate2023,
    remap = lookupSitenames2023[, c(2,10,6)], # remap needs to be old, new, type
    repeatSamples = TRUE,
    fillColor = c(
        supportingColorPalette[c(1, 5, 2, 6, 7, 3, 8)],
        corporateColorPalette[c(1, 2)],
        supportingColorPalette[c(4,1,2,3,4,5)]
    ),
    levelsPlot = levelOrder
) + theme(legend.position = "bottom")





## For each site (pooled replicates for each land use and site)

aitchisonSite2023 <- vegdist(compSite2023, "euc", diag = F)

min(aitchisonSite2023)
max(aitchisonSite2023)

## create plots

siteHeatmap <-
    aitchisonSite2023 %>% aitHeatmap(fillColor1 = supportingColorPalette[2],
                                      fillColor2 = corporateColorPalette[4]) 

lookupSitenames2023$treatment.site <- paste0(lookupSitenames2023$treatment.code.original, "-S", lookupSitenames2023$site.code)

levelOrder = c(
    "CF-CF",
    "CF-CF2",
    "CF-I",
    "CF-IA",
    "CF-IB",
    "CF-R",
    "CF2-CF2",
    "CF2-I",
    "CF2-IA",
    "CF2-IB",
    "CF2-R",
    "I-I",
    "I-IA",
    "I-IB",
    "I-R",
    "IA-IB",
    "IA-R",
    "IB-R",
    "R-R"
)

aitComparison(
    inputDist = aitchisonSite2023,
    remap = lookupSitenames2023[, c(11,11,6)], # remap needs to be old, new, type
    repeatSamples = TRUE,
    fillColor = c(
        supportingColorPalette[c(1, 5, 2, 6, 7, 3, 8)],
        corporateColorPalette[c(1, 2, 3:6)],
        supportingColorPalette
    ),
    levelsPlot = levelOrder
) + theme(legend.position = "bottom")



## For each land cover type (pooled replicates across sites)

aitchisonLC2023 <- vegdist(compType2023, "euc", diag = F)

min(aitchisonLC2023)
max(aitchisonLC2023)

## create plots

typeHeatmap <-
    aitchisonLC2023 %>% aitHeatmap(fillColor1 = supportingColorPalette[2],
                                     fillColor2 = corporateColorPalette[4]) 


ggsave("RTerraDistHeatmap_2023.pdf",
       plot = typeHeatmap,
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300
)

remove(siteHeatmap, typeHeatmap)


## ----- Indicator 3: PCA qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/


# create replicate pca
pca_repl <- compMatrix2023 %>%
    PCA(., scale.unit = F, graph = F)

fviz_pca_ind(
    pca_repl,
    geom.ind = "point",
    col.ind = str_split_fixed(rownames(compMatrix2023), "-",5)[,4],
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "LC Type",
    palette = supportingColorPalette[c(1:8)]
)


# create LC-Site pca
pca_LCplots <- compSite2023 %>%
    PCA(., scale.unit = F, graph = F)

groupstemp <- str_split_fixed(rownames(compSite2023), "-",2)[,1]
groupstemp[groupstemp %in% c("IA", "IB")] <- "I"

viz_pcaLCPlots <- fviz_pca_ind(
    pca_LCplots,
    geom.ind = "point",
    col.ind = groupstemp,
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "LC Type",
    palette = supportingColorPalette[c(1:8)]
)

ggpubr::ggpar(viz_pcaLCPlots,
              title = "Community Composition Visualization using PCA",
              subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F
)

ggsave("RTerraPCA_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)


viz_pcaPlots_contrib <- fviz_contrib(pca_LCplots, choice = "ind", axes = 1:2)

# fviz_pca_biplot(pca_LCplots, 
#                 # Sites 
#                 col.ind = groupstemp,
#                 addEllipses = T,
#                 ellipse.type = "convex",
#                 label = "var",
#                 repel = T,
#                 alpha.var ="contrib")

remove(pca_plots, viz_pcaPlots, viz_pcaPlots_contrib)
