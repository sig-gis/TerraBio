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
length(apuiASV$ASVHeader[apuiASV$sample == "EM110-01A-R1" &
                             apuiASV$readOrigin == "merged" &
                             apuiASV$relativeAbundanceOnSample > 0.05])


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


# Create accumulation curve for all replicates or sites sampled?


## ----- Indicator 1: Alpha Diversity -------------------------------------------
# I1: Alpha diversity using Effective Species Richness

# Create a table with the alpha diversity measures for each replicate. 
sampleNumbers <- str_split(rownames(apuiMatrix), pattern = "-", simplify = T)[ , 3] %>%
    substr(.,1,2) %>%     factor(., levels = unique(.))

sampleNames <- str_split(rownames(apuiMatrix), pattern = "-", simplify = T)[ , c(1,2)]
    sampleNames <- paste0(sampleNames[,1], "-", sampleNames[,2]) %>% factor(., levels = unique(.))

sampleTreatments <- str_split(rownames(apuiMatrix), pattern = "-", simplify = T)[ , 2] %>%
    factor(., levels = unique(.))

## combine replicates for each sample, then calculate alpha for each sample. can
## use "group" for this.

# temp <- specnumber(apuiMatrix, MARGIN = 1)
# plot(temp,
#      pch = 19,
#      col = factor(sampleTreatments))
# legend(x = 40, y = 300,
#        legend = c("Forest", "Succession", "Pasture", "SAF"),
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
#        legend = c("Forest", "Succession", "Pasture", "SAF"),
#        fill = c("black", "green", "red", "blue"),
#        cex = 1)
# 
# plot(temp2,
#      pch = 19,
#      col = factor(sampleTreatments[seq(1, length(sampleTreatments), 3)]))
# legend(x = 12, y = 510,
#        legend = c("Forest", "Succession", "Pasture", "SAF"),
#        fill = c("black", "green", "red", "blue"),
#        cex = 1)

# Some make sense, but in general the patterns are very off. for example, Succession around 250 for each replicate, but then only 110 for group function.
# Not sure that this is my issue. I think this is with the new group functionality for vegan... can try updating R, RStudio and packages and then rerunning code.

# This was a problem with the group function, namely issues caused because aggregate reorders, silently changing the output order. I changed my summary function to fix it.




# Compare diversities between land use types, one dot for each farm field sampled


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
    distinct() %>%
    arrange(siteType) %>%
    mutate(
        treatmentYear = c("*", "", "", "2022/3", "", "2021", "2022/3", "", "2019", "", "2022/3", "2019", "", "2014", "", "")
    )




# Plot for ESR

apuiSampleAlpha %>%
    mutate(treatment = factor(treatment, levels = c("Pasture", "SAF", "Forest", "Succession"))) %>%
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

apuiSampleAlpha %>%
    mutate(treatment = factor(treatment, levels = c("Pasture", "SAF", "Forest", "Succession"))) %>%
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
    scale_color_manual(values = supportingColorPalette) +
    geom_text(aes(label = treatmentYear), nudge_x = .2)

ggsave("ApuiDateESR_2022.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
       )


sampleAlphaSummary <- apuiSampleAlpha %>% 
    group_by(treatment) %>%
    summarise(n = n(),
              meanESR = mean(effectiveSR))

# Compare diversities between land use types, pooled for all fields of a land
# use type sampled--note this isn't great because of the different number of
# fields sampled...

apuiTreatmentAlpha <- alphaGroupMetrics(apuiMatrix,
                                        groupNames = sampleTreatments) 

apuiTreatmentAlpha %>%
    mutate(siteType = factor(
        siteType,
        levels = c("Pasture", "SAF", "Forest", "Succession")
    )) %>%
    ggplot(aes(siteType, effectiveSR, fill = siteType)) +
    geom_bar(stat = "identity") +
    labs(fill = "Site Type",
         y = "Effective Species Richness",
         x = "Site Type") +
    theme(legend.position = "bottom"
          ) +
    scale_fill_manual(values = supportingColorPalette) +
    geom_text(stat = 'identity',
              aes(label = round(effectiveSR)),
              position = position_stack(vjust = 0.5))




## ----- Convert matrices to compositional data ------------------

# create compositional matrices
# This is for each replicate (no pooling)
apuiCompMatrix <- compMatrix(apuiMatrix)

# This is for pooled farm fields (samples)
apuiSampleCompMatrix <- apuiMatrix %>%
    aggregate.data.frame(., list(sampleNames), sum)
    rownames(apuiSampleCompMatrix) <- apuiSampleCompMatrix$Group.1 
    apuiSampleCompMatrix <- compMatrix(apuiSampleCompMatrix[ , -1])

# This is for pooled treatments
apuiTreatmentCompMatrix <- apuiMatrix %>%
    aggregate.data.frame(., list(sampleTreatments), sum)
    rownames(apuiTreatmentCompMatrix) <- apuiTreatmentCompMatrix$Group.1 
    apuiTreatmentCompMatrix <- compMatrix(apuiTreatmentCompMatrix[ , -1])
    

    
    
    
apuiPairedCompMatrix <- compMatrix(apuiPairedMatrix) # this is for each replicate, to be used for PI 4



## ----- Indicator 2: Beta diversity w/ Aitchison distance -----------------------------
# Aitchison distance uses the Euclidian distance of the compositional data that
# has been center log transformed; see
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/

# For each replicate--internal use only, really.
aitchisonReplicate <- vegdist(apuiCompMatrix, "euc", diag = F)

min(aitchisonReplicate)
max(aitchisonReplicate)

aitchisonReplicate %>% aitHeatmap()


## For each sample (pooled replicates by farm field)
aitchisonSample <- vegdist(apuiSampleCompMatrix, "euc", diag = F)

min(aitchisonSample)
max(aitchisonSample)

aitchisonSample %>% aitHeatmap() # note this will have some in strange positions--it's fine, we're not using this for anything other than diagnostics.

levelOrder = c("Pasture-Pasture", "Forest-Forest", "Succession-Succession", "SAF-SAF",
               "Pasture-Succession", "Pasture-SAF",
               "Forest-Succession", "Forest-SAF", "Forest-Pasture", 
               "SAF-Succession")


remap <- lookupSitenames %>%
    transmute(
        shortSampleID = shortSampleID,
        prettyName = shortSampleID,
        Treatment = treatment
    ) %>%
    distinct()

aitComparison(
    inputDist = aitchisonSample,
    remap = remap, 
    levelsPlot = levelOrder,
    repeatSamples = TRUE,
    fillColor = c(
        supportingColorPalette[c(3, 5, 6, 7, 1)],
        corporateColorPalette[c(2:3)],
        supportingColorPalette[c(2,8,4)]
        
    )
)

ggsave("ApuiDistBoxplots_2022.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)

## ----- PI3: Change in beta diversity -----------------------------------------

# Looking at the overall distance between the treatments. 


#aggregate by site 




aitchisonTreatment <- vegdist(apuiTreatmentCompMatrix[ order(rownames(apuiTreatmentCompMatrix)), ], "euc") 

min(aitchisonTreatment)
max(aitchisonTreatment)

## create plots

treatmentHeatmap <-
    aitchisonTreatment %>% aitHeatmap(fillColor1 = supportingColorPalette[2],
                                      fillColor2 = corporateColorPalette[4]) 


ggsave("ApuiDistHeatmap_2022.pdf",
       plot = treatmentHeatmap,
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 6,
       units = "in",
       dpi = 300
)

## ----- Indicator 4: qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# create replicate pcas
pca_replicates <- apuiCompMatrix %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaReplicates <- fviz_pca_ind(
    pca_replicates,
    geom.ind = "point",
    col.ind = str_split(rownames(apuiCompMatrix), pattern = "-", simplify = T)[ , 2],
    addEllipses = T,
    ellipse.type = "confidence",
    pointsize = 3,
    mean.point = F,
    legend.title = "Site Type",
    palette = supportingColorPalette[c(3,1,2,4)]
)
ggpubr::ggpar(viz_pcaReplicates,
              title = "Community Composition Visualization using PCA -- Replicates",
              subtitle = paste0(phylum, collapse = " ")#, xlab = F, ylab = F, tickslab = F
)

viz_pcaPlots_repli <- fviz_contrib(pca_replicates, choice = "ind", axes = 1:2)


# create field sample pcas (pooled replicates)
pca_samples <- apuiSampleCompMatrix %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaSamples <- fviz_pca_ind(
    pca_samples,
    geom.ind = "point",
    col.ind = str_split(rownames(apuiSampleCompMatrix), pattern = "-", simplify = T)[ , 2],
    addEllipses = T,
    repel = T,
    ellipse.type = "confidence",
    ellipse.level = 0.95,
    pointsize = 3,
    mean.point = F,
    legend.title = "Site Type",
    palette = supportingColorPalette[c(3,1,2,4)], axes = c(1,2)
)
ggpubr::ggpar(viz_pcaSamples,
              title = paste0("Community Composition Visualization using PCA"),#, hortaSubset),
              subtitle = paste0(phylum, collapse = " "),
              xlab = F, ylab = F, tickslab = F
)

viz_pcaPlots_contrib <- fviz_contrib(pca_plots, choice = "ind", axes = 1:2)



ggsave("ApuiPCA_2022.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "OutputImages/",
       width = 8,
       height = 5,
       units = "in",
       dpi = 300
)
