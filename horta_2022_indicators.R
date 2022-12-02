# This analysis is based on soil samples collected in June at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and determines the
# TerraBio biodiversity indicators. There are also some species accumulation
# curves to verify sampling completeness.

## ----- Data ingestion & setup -----------------------------------

# libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(iNEXT)

source("allianceBranding.R")
source("functions.R")

source("../../../RCode/R_Scripts/triplet_fixer.R") # From my github repository
source("../../../RCode/R_Scripts/repeat_multipatt.R") # ditto

source("horta_2022_data_processing.R")
remove(list = grep("volume",ls(), value = T))
remove(abundanceLetter)

rare = 50
hortaSubset = "b-R1" # indicates we're using the buffer/silica samples for R1 primer set.
hortaGroups <- c(rep("Counterfactual", 3), rep("Forest", 3), rep("Restoration", 3), rep("Syntropic", 3))



hortaMatrixSmall <-
    hortaMatrix[grepl(pattern = hortaSubset, rownames(hortaMatrix)),]
hortaMatrixSmall <-
    hortaMatrixSmall[, colSums(hortaMatrixSmall) > 0]

hortaLetterSmall <-
    hortaMatrixLetter[grepl(pattern = hortaSubset, rownames(hortaMatrixLetter)),]
hortaLetterSmall <-
    hortaLetterSmall[, colSums(hortaLetterSmall) > 0]

hortaMatrixR1 <-
    hortaMatrix[grepl(pattern = "R1", rownames(hortaMatrix)),]
hortaMatrixR1 <-
    hortaMatrixR1[, colSums(hortaMatrixR1) > 0]


hortaCompLetter <- compMatrix(inputMatrix = hortaLetterSmall)

hortaCompMatrix <- compMatrix(inputMatrix = hortaMatrixSmall)

hortaCompR1 <- compMatrix(inputMatrix = hortaMatrixR1)

## ----- Storage method graphs --------------------------------

# plot amount of DNA recovered from each sample by buffer/silica
ggplot(hortaInfo, aes(x = storage, y = concentrationDNA_nguL)) + 
    geom_boxplot() + facet_grid(. ~ volumeSampleID)

# plot 'quality' measure for each sample by buffer/silica
ggplot(hortaInfo, aes(x = storage, y = purityDNA)) + 
    geom_boxplot() + facet_grid(. ~ volumeSampleID)

# plot all ASV absolute abundance by buffer/silica
ggplot(hortaASV, aes(x = preservation, y = asvAbsoluteAbundance)) + 
    geom_boxplot() + geom_jitter()

# plot total ASV absolute abundance grouped by sample for buffer/silica
hortaSampleASV <- hortaASV %>% group_by(metadata_4, primerName, preservation) %>% 
    summarise(totalAbund = sum(asvAbsoluteAbundance), countASV = length(unique(ASVHeader)))

ggplot(hortaSampleASV, aes(x = preservation, y = totalAbund)) + 
    geom_boxplot() + geom_jitter(aes(color = primerName))

ggplot(hortaSampleASV, aes(x = preservation, y = totalAbund)) + 
    geom_boxplot() + geom_jitter() + 
    facet_grid(primerName ~ ., scales = "free")

# plot total number of ATVs grouped by sample for buffer/silica
ggplot(hortaSampleASV, aes(x = preservation, y = countASV)) + 
    geom_boxplot() + geom_jitter(aes(color = primerName))

ggplot(hortaSampleASV, aes(x = preservation, y = countASV)) + 
    geom_boxplot() + geom_jitter() + 
    facet_grid(primerName ~ ., scales = "free")


## ----- Accumulation curves ----------------------------------

# create species accumulation curves for both buffer and silica (pretending 6
# replicates instead of 3)
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-CF1|R1-CF2|R1-CF3",
                                   x = rownames(hortaMatrix)), ]),
     xlab = "Number of replicates [R1]",
     ylab = "Number of species", ylim = c(10,450))
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-Co1|R1-Co2|R1-Co3",
                                   x = rownames(hortaMatrix)), ]),
     add = TRUE, col = "green")
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-Re1|R1-Re2|R1-Re3",
                                   x = rownames(hortaMatrix)), ]),
     add = TRUE, col = "blue")
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-Sy1|R1-Sy2|R1-Sy3",
                                   x = rownames(hortaMatrix)), ]),
     add = TRUE, col = "red")
legend(x = 1, y = 400,
       legend = c("Counterfactual", "Control", "Restoration", "Syntropic"),
       fill = c("black", "green", "blue", "red"),
       cex = 1)


# create species accumulation curves for both buffer and silica (pretending 6
# replicates instead of 3); remove rare species
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-CF1|R1-CF2|R1-CF3",
                                   x = rownames(hortaMatrix)), 
                             colSums(hortaMatrix) > rare]),
     xlab = "Number of replicates [R1, rare = 50]",
     ylab = "Number of species", ylim = c(10,200))
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-Co1|R1-Co2|R1-Co3",
                                  x = rownames(hortaMatrix)),
                            colSums(hortaMatrix) > rare]),
     add = TRUE, col = "green")
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-Re1|R1-Re2|R1-Re3",
                                  x = rownames(hortaMatrix)),
                            colSums(hortaMatrix) > rare]),
     add = TRUE, col = "blue")
plot(specaccum(hortaMatrix[ grepl(pattern = "R1-Sy1|R1-Sy2|R1-Sy3",
                                  x = rownames(hortaMatrix)),
                            colSums(hortaMatrix) > rare]),
     add = TRUE, col = "red")
legend(x = 1, y = 200,
       legend = c("Counterfactual", "Control", "Restoration", "Syntropic"),
       fill = c("black", "green", "blue", "red"),
       cex = 1)



## ----- Proposed Indicator 1: Alpha -------------------------------------------
# PI1: Alpha diversity

# Create a table with the alpha diversity measures for each replicate. 
## *** NEED TO ADD FIELD (REPLICATE)***
hortaAlpha<- alphaMetrics(hortaMatrixSmall, 
                          groupNames = hortaGroups, 
                          replNames = str_sub(rownames(hortaMatrixSmall), -3,-1))
hortaAllAlpha<- alphaMetrics(hortaMatrixR1,
                             groupNames = rep(hortaGroups, 2),
                             replNames = str_sub(rownames(hortaMatrixR1), -3,-1))

hortaAlphaGroup <- alphaGroupMetrics(hortaMatrixSmall,
                                     groupNames = hortaGroups)


# Test to compare site diversities between land use types

# Test for raw species number
# 
# speciesRLMER <- lme4::lmer( speciesRichness ~
#                                   siteType +
#                                   (1 | siteType),
#                               data = hortaAlpha,
#                               REML = TRUE)
# anova(speciesRLMER)
# lmerTest::rand(speciesRLMER)
# summary(speciesRLMER)
# test_speciesRLMER<-car::Anova(mod = speciesRLMER)
# emmeans::emmeans(speciesRLMER, pairwise~siteType)    
#   
# 
# hortaAlpha %>%
#     ggplot() +
#     geom_boxplot(aes(siteType, speciesRichness),
#                  outlier.shape = NA) +
#     geom_jitter(
#         aes(siteType, speciesRichness, color = siteType),
#         width = 0.1,
#         height = 0
#     ) +
#     labs(color = "Site Type", y = "Raw Species Richness",
#          x = "Site Type") + 
#     theme(legend.position="bottom")


# Test for ESR


hortaAlpha %>%
    mutate(siteType = factor(siteType, levels = c("Counterfactual", "Syntropic", "Forest", "Restoration"))) %>%
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
    theme(legend.position="bottom",
          text = element_text(family = "Calibri"))+
    scale_color_manual(values = supportingColorPalette)

hortaAlphaGroup %>% 
    mutate(siteType = factor(siteType, levels = c("Counterfactual", "Syntropic", "Forest", "Restoration"))) %>%
    ggplot(
        aes(siteType, effectiveSR, fill = siteType)
        ) +
    geom_bar(stat = "identity") +
    labs(fill = "Site Type", 
         y = "Effective Species Richness",
         x = "Site Type") +
    theme(legend.position = "bottom",
          text = element_text(family = "Calibri")) +
    scale_fill_manual(values = supportingColorPalette) +
    geom_text(stat='identity', aes(label=round(effectiveSR)),position = position_stack(vjust = 0.5))
    


# # Test for inverse Simpson
# 
# invSimpsonLMER <- lme4::lmer(invSimpson ~
#                                  siteType +
#                                  (1 | siteType),
#                              data = hortaAlpha,
#                              REML = TRUE)
# anova(invSimpsonLMER)
# lmerTest::rand(invSimpsonLMER)
# summary(invSimpsonLMER)
# test_shannonLMER <- car::Anova(invSimpsonLMER)
# emmeans::emmeans(invSimpsonLMER, pairwise~siteType)    
# 
# hortaAlpha %>%
#     ggplot() +
#     geom_boxplot(aes(siteType, invSimpson),
#                  outlier.shape = NA) +
#     geom_jitter(
#         aes(siteType, invSimpson, color = siteType),
#         width = 0.1,
#         height = 0
#     ) +
#     labs(color = "Field Type", y = "Inv. Simpson Diversity",
#          x = "Site Type") +
#     theme(legend.position="bottom")


## ----- PI2: Beta diversity w/ Aitchison distance -----------------------------
# aitchison distance uses the euclidian distance of the compositional data that
# has been center log transformed; see
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/


aitchisonPlot <- vegdist(hortaCompMatrix, "euc", diag = F)

min(aitchisonPlot)
max(aitchisonPlot)



## create plots

#aitchisonPlot %>% aitHeatmap()

levelOrder = c("Counterfactual-Counterfactual", "Forest-Forest", "Restoration-Restoration", "Syntropic-Syntropic",
               "Forest-Counterfactual", "Restoration-Counterfactual", "Syntropic-Counterfactual",
               "Restoration-Forest", "Syntropic-Forest",
               "Syntropic-Restoration")
aitComparison(aitchisonPlot, remap = lookupSitenames, repeatSamples = TRUE, fillColor = c(supportingColorPalette[c(1,5,2,6,7,3,8)], corporateColorPalette[c(1,2)], supportingColorPalette[4]), levelsPlot = levelOrder)


    

## ----- PI3: Change in beta diversity -----------------------------------------

# Looking at the overall distance between the treatments. 


#aggregate by site 

aitchisonSite <- vegdist(hortaCompLetter, "euc") 

min(aitchisonSite)
max(aitchisonSite)

## create plots

siteHeatmap <- aitchisonSite %>% aitHeatmap(fillColor1 = supportingColorPalette[2], fillColor2 = corporateColorPalette[4])
siteHeatmap +
    scale_x_discrete(labels=c("b-R1-Re" = "Restoration",
                              "b-R1-Co" = "Forest",
                              "b-R1-Sy" = "Syntropic",
                              "b-R1-CF" = "Counterfactual",
                              "s-R1-Re" = "Restoration",
                              "s-R1-Co" = "Forest",
                              "s-R1-Sy" = "Syntropic",
                              "s-R1-CF" = "Counterfactual")) +
    scale_y_discrete(labels=c("b-R1-Re" = "Restoration",
                              "b-R1-Co" = "Forest",
                              "b-R1-Sy" = "Syntropic",
                              "b-R1-CF" = "Counterfactual",
                              "s-R1-Re" = "Restoration",
                              "s-R1-Co" = "Forest",
                              "s-R1-Sy" = "Syntropic",
                              "s-R1-CF" = "Counterfactual"))

#aitComparison(aitchisonSite, lookupSitenames, FALSE)


## ----- PI 4: qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# create plot pcas
pca_plots <- hortaCompMatrix %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaPlots <- fviz_pca_ind(
    pca_plots,
    geom.ind = "point",
    col.ind = hortaGroups,
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "Site Type",
    palette = supportingColorPalette[c(1,3,4,2)]
)
ggpubr::ggpar(viz_pcaPlots,
              title = paste0("Community Composition Visualization using PCA"),#, hortaSubset),
              subtitle = paste0(phylum, collapse = " "), xlab = F, ylab = F, tickslab = F
              )

viz_pcaPlots_contrib <- fviz_contrib(pca_plots, choice = "ind", axes = 1:2)

fviz_pca_biplot(pca_plots,
                # Sites
                col.ind = hortaGroups,
                addEllipses = T,
                ellipse.type = "convex",
                label = "var",
                repel = T,
                max.overlaps = 5,
                alpha.var ="contrib")

# create plot pcas for both 
# pca_plotsAll <- hortaCompR1 %>%
#     PCA(., scale.unit = F, graph = F)
# 
# viz_pcaPlotsAll <- fviz_pca_ind(
#     pca_plotsAll,
#     geom.ind = "point",
#     col.ind = c(paste0(hortaGroups, "-b"), paste0(hortaGroups, "-s")),
#     addEllipses = T,
#     ellipse.type = "convex",
#     legend.title = "Group"
# )
# ggpubr::ggpar(viz_pcaPlotsAll,
#               title = paste0("Plots - PCA - ", paste0(phylum, collapse = " "), " - R1"))

