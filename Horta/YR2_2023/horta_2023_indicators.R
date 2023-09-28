# This analysis is based on soil samples collected at Horta and processed by
# EcoMol. There are four land uses: CF = counter factual, Co = control/forest,
# Re = restauracao, Sy = syntropic.

# This analysis includes data from 2022 and 2023

# This file take the eDNA species table provided by EcoMol and determines the
# TerraBio biodiversity indicators. There are also some species accumulation
# curves to verify sampling completeness.

# Code written Sept. 2023 by Karen Dyson

# to do list:
#    + Four indicators--which need special multi-year functions?




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

source("horta_2023_data_processing.R")
remove(abundanceLetter)

rare = 50
hortaGroups2022 <- c(rep("Counterfactual", 3), rep("Forest", 2),
                     rep("Restoration", 3), rep("Syntropic", 3))
hortaGroups2023 <- c(rep("Forest", 3), rep("Restoration", 3),
                     rep("Syntropic", 3), rep("Counterfactual", 3))



# Create compositional matrices

compMatrix2022 <- compMatrix(inputMatrix = hortaMatrix2022)

compMatrix2023 <- compMatrix(inputMatrix = hortaMatrix2023)

compType2022 <- compMatrix(inputMatrix = hortaMatrixType2022)

compType2023 <- compMatrix(inputMatrix = hortaMatrixType2023)

# the 2023 matrices had adjusted imputations--check


## ----- Proposed Indicator 1: Alpha -------------------------------------------
# PI1: Alpha diversity

# Create a table with the alpha diversity measures for each replicate. 



hortaAlpha2022 <- alphaMetrics(hortaMatrix2022, 
                          groupNames = hortaGroups2022, 
                          replNames = str_sub(rownames(hortaMatrix2022), -3,-1))
hortaAlpha2023 <- alphaMetrics(hortaMatrix2023,
                             groupNames = hortaGroups2023,
                             replNames = str_sub(rownames(hortaMatrix2023), -2,-1))

hortaAlphaGroup2022 <- alphaGroupMetrics(hortaMatrixType2022,
                                     groupNames = c("Counterfactual", "Forest", "Restoration", "Syntropic"))
hortaAlphaGroup2023 <- alphaGroupMetrics(hortaMatrixType2023,
                                         groupNames = c("Forest", "Restoration", "Syntropic", "Counterfactual"))


# ESR for restoration in 2023 is weirdly low.


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
    theme(legend.position = "bottom") +
    scale_fill_manual(values = supportingColorPalette) +
    geom_text(stat='identity', aes(label=round(effectiveSR)),position = position_stack(vjust = 0.5))

ggsave("HortaESR_2022.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)
