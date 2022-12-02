# This analysis is based on soil samples collected in September 2022 at Cafe Apui
# and processed by EcoMol.

# This file take the eDNA species table provided by EcoMol and determines the
# TerraBio biodiversity indicators. There are also some species accumulation
# curves to verify sampling completeness.

## ----- Data ingestion & setup -----------------------------------

# libraries
library(iNEXT)

source("apui_2022_data_processing.R")

# checks
length(apuiASV$ASVHeader[apuiASV$sample == "EM110-01A-R1" & apuiASV$readOrigin == "merged" & apuiASV$relativeAbundanceOnSample > 0.05])


## ----- Convert matrices to compositional data ------------------

# remove columns with 0s that might result from split
apuiMatrix <- apuiMatrix [ , colSums(apuiMatrix) > 0 ]
apuiPairedMatrix <- apuiPairedMatrix[ , colSums(apuiPairedMatrix) > 0 ]

# create compositional matrices
apuiCompMatrix <- compMatrix(apuiMatrix)
apuiPairedCompMatrix <- compMatrix(apuiPairedMatrix)



## ----- Accumulation curves ----------------------------------

# create species accumulation curves for both buffer and silica for farm 3
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-Pasture",
                                  x = rownames(apuiPairedMatrix)), ]),
     xlab = "Number of replicates at Farm 3, both Buffer and Silica",
     ylab = "Number of species",
     ylim = c(0,800))
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-Forest",
                                  x = rownames(apuiPairedMatrix)), ]),
     add = TRUE, col = "green")
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-SAF",
                                  x = rownames(apuiPairedMatrix)), ]),
     add = TRUE, col = "blue")

legend(x = 1, y = 600,
       legend = c("Pasture", "Forest", "SAF"),
       fill = c("black", "green", "blue"),
       cex = 1)


# create species accumulation curves for each buffer and silica
#Buffer
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-Pasture",
                                       x = rownames(apuiPairedMatrix)) &
                                     grepl(pattern = "Buffer",
                                           x = rownames(apuiPairedMatrix)), ]),
     xlab = "Number of replicates at Farm 3, split Buffer and Silica",
     ylab = "Number of species",
     ylim = c(0,550))
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-Forest",
                                       x = rownames(apuiPairedMatrix)) &
                                     grepl(pattern = "Buffer",
                                           x = rownames(apuiPairedMatrix)), ]),
     add = TRUE, col = "green")
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-SAF",
                                       x = rownames(apuiPairedMatrix)) &
                                     grepl(pattern = "Buffer",
                                           x = rownames(apuiPairedMatrix)), ]),
     add = TRUE, col = "blue")

# Silica
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-Pasture",
                                       x = rownames(apuiPairedMatrix)) &
                                     grepl(pattern = "Silica",
                                           x = rownames(apuiPairedMatrix)), ]),
     add = TRUE, lty = 2)
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-Forest",
                                       x = rownames(apuiPairedMatrix)) &
                                     grepl(pattern = "Silica",
                                           x = rownames(apuiPairedMatrix)), ]),
     add = TRUE, col = "green", lty = 2)
plot(specaccum(apuiPairedMatrix[ grepl(pattern = "Farm3-SAF",
                                       x = rownames(apuiPairedMatrix)) &
                                     grepl(pattern = "Silica",
                                           x = rownames(apuiPairedMatrix)), ]),
     add = TRUE, col = "blue", lty = 2)


legend(x = 1, y = 500,
       legend = c("Pasture", "Forest", "SAF"),
       fill = c("black", "green", "blue"),
       cex = 1)


## ----- Proposed Indicator 1: Alpha -------------------------------------------
# PI1: Alpha diversity

# Create a table with the alpha diversity measures for each replicate. 
sampleNumbers <- str_split(rownames(apuiMatrix), pattern = "-", simplify = T)[ , 3] %>%
    substr(.,1,2)
sampleTreatments <- str_split(rownames(apuiMatrix), pattern = "-", simplify = T)[ , 2]

## combine replicates for each sample, then calculate alpha for each sample. can
## use "group" for this I believe.




# temp <- specnumber(apuiMatrix, MARGIN = 1)
# plot(temp,
#      pch = 19,
#      col = factor(sampleTreatments))
# legend(x = 40, y = 300,
#        legend = c("Forest", "Regeneration", "Pasture", "SAF"),
#        fill = c("black", "green", "red", "blue"),
#        cex = 1)
# 
# temp <- specnumberMOD(x = apuiMatrix, MARGIN = 1, groups = sampleNumbers)
# temp2 <- specnumber(x = apuiMatrix, MARGIN = 1, groups = sampleNumbers)
# 
# plot(temp,
#      pch = 19,
#      col = factor(sampleTreatments[seq(1, length(sampleTreatments), 3)]))
# legend(x = 12, y = 510,
#        legend = c("Forest", "Regeneration", "Pasture", "SAF"),
#        fill = c("black", "green", "red", "blue"),
#        cex = 1)
# 
# plot(temp2,
#      pch = 19,
#      col = factor(sampleTreatments[seq(1, length(sampleTreatments), 3)]))
# legend(x = 12, y = 510,
#        legend = c("Forest", "Regeneration", "Pasture", "SAF"),
#        fill = c("black", "green", "red", "blue"),
#        cex = 1)

# Some make sense, but in general the patterns are very off. for example, regeneration around 250 for each replicate, but then only 110 for group function.
# Not sure that this is my issue. I think this is with the new group functionality for vegan... can try updating R, RStudio and packages and then rerunning code.

# If it still fails after these measures, I'll need to pool replicate data for the samples and then do the calculations. 


# Plot for raw species number

temp %>%
    ggplot() +
    geom_boxplot(aes(treatment, speciesRichness),
                 outlier.shape = NA) +
    geom_jitter(
        aes(treatment, speciesRichness, color = treatment),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Site Type", y = "Raw Species Richness",
         x = "Site Type") + 
    theme(legend.position="bottom")







apuiSampleAlpha<- alphaGroupMetrics(apuiMatrix, 
                          groupNames = sampleNumbers) %>%
    mutate(
        siteType = paste0("EM110-", siteType)
    ) %>%
    left_join(
        .,
        y = lookupSitenames[ , colnames(lookupSitenames) %in% c("farmPerson", "farmName", "treatment", "EM.ID")],
        by = c("siteType"="EM.ID")
    ) %>%
    distinct()



# Plot to compare site diversities between land use types

# Plot for ESR

apuiSampleAlpha %>%
    mutate(treatment = factor(treatment, levels = c("Pasture", "SAF", "Forest", "Regeneration"))) %>%
    ggplot(aes(treatment, effectiveSR)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = treatment),
        width = 0.05,
        height = 0,
        size = 4
    ) +
    labs(color = "Site Type", y = "Effective Species Richness",
         x = "Site Type") +
    theme(legend.position = "bottom"
          ) +
    scale_color_manual(values = supportingColorPalette) 



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



