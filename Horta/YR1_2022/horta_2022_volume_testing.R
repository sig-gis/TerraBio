# This analysis is based on soil samples collected in June at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and determines three
# things: 
# 1. Is Longmire's buffer or desiccant a better preservative to use in
# the field, where "better" is more community data and coverage returned?
# 2. How many replicates are needed per sample? 
# 3. How much soil should be collected from each site?

## ----- Data ingestion ---------------------------------------

# libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(iNEXT)

source("functions.R")
source("../../../RCode/R_Scripts/triplet_fixer.R")
source("horta_2022_data_processing.R")

rare = 50


## ----- Vol Testing | Preservative --------------------------
# plot all ASV absolute abundance by buffer/silica
ggplot(volumeASV, aes(x = preservation, y = asvAbsoluteAbundance)) + 
    geom_boxplot() + geom_jitter()

# plot total ASV absolute abundance grouped by sample for buffer/silica
volumeSampleASV <- volumeASV %>% group_by(metadata_4, primerName, preservation) %>% 
    summarise(totalAbund = sum(asvAbsoluteAbundance), countASV = length(unique(OTU)))

ggplot(volumeSampleASV, aes(x = preservation, y = totalAbund)) + 
    geom_boxplot() + geom_jitter(aes(color = primerName))

ggplot(volumeSampleASV, aes(x = preservation, y = totalAbund)) + 
    geom_boxplot() + geom_jitter() + 
    facet_grid(primerName ~ ., scales = "free")

# plot total number of OTUs grouped by sample for buffer/silica
ggplot(volumeSampleASV, aes(x = preservation, y = countASV)) + 
    geom_boxplot() + geom_jitter(aes(color = primerName))

ggplot(volumeSampleASV, aes(x = preservation, y = countASV)) + 
    geom_boxplot() + geom_jitter() + 
    facet_grid(primerName ~ ., scales = "free")


# plot amount of DNA recovered from each sample by buffer/silica
ggplot(volumeInfo, aes(x = storage, y = concentrationDNA_nguL)) + 
    geom_boxplot()# + geom_jitter(aes(color = volumeSampleID))

# plot 'quality' measure for each sample by buffer/silica
ggplot(volumeInfo, aes(x = storage, y = purityDNA)) + 
    geom_boxplot() #+ geom_jitter(aes(color = volumeSampleID))

## ----- Bonus | Primers ---------------------------------------
# What's the overlap between the communities for buffer?
table(volumeASV$order_BLASTn[volumeASV$primerName == "R1"])
table(volumeASV$genusBLASTn[volumeASV$primerName == "R2"])

# plot the relative number of each order for R1 and R2
ggplot(volumeASV, aes(x = order_BLASTn)) + 
    geom_bar(aes(fill = primerName), position = "fill") +
    coord_flip()

genusR2 <- volumeASV$genusBLASTn[volumeASV$primerName == "R2"]
genusR1 <- volumeASV$genusBLASTn[volumeASV$primerName == "R1"]

genusR2[!(genusR2 %in% genusR1)]

remove(genusR1, genusR2)


#Create a compositional matrix

# square-root Bayesian-multiplicative replacement of zeros with the cmultRepl()
# function (Ladin et al., 2021) 

# Most other code in GitHub uses CZM. e.g. See
# https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-Biplot
# and https://raw.githubusercontent.com/ggloor/CoDaSeq/6ff864aade46cd3c8b0eff3bb54d5460775f92cd/CoDaSeq/vignettes/CoDaSeq_vignette.Rnw
# This latter contends that this is the most principled method.
# Also see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/

# For the buffer
volumeMatrixB <- volumeMatrix[grepl(pattern = "b-R",
                                   x = rownames(volumeMatrix)),]
volumeMatrixB <- volumeMatrixB[ , colSums(volumeMatrixB) > 0]

volumeMatrixB_0repl <-
    cmultRepl(volumeMatrixB,
              label = 0,
              method = "CZM")

volumeMatrixB_comps <- cdt.acomp(x = volumeMatrixB_0repl) %>% 
    as_tibble(., rownames = NA)

# For the silica
volumeMatrixS <- volumeMatrix[grepl(pattern = "s-R",
                                    x = rownames(volumeMatrix)),]
volumeMatrixS <- volumeMatrixS[ , colSums(volumeMatrixS) > 0]

volumeMatrixS_0repl <-
    cmultRepl(volumeMatrixS,
              label = 0,
              method = "CZM")

volumeMatrixS_comps <- cdt.acomp(x = volumeMatrixS_0repl) %>% 
    as_tibble(., rownames = NA)

#cleanup
remove(volumeMatrixB_0repl, volumeMatrixS_0repl)

# Now calculate Aitchison distances.
aitchisonPlotB <- vegdist(volumeMatrixB_comps, method = "eucl", diag = F)
aitchisonPlotS <- vegdist(volumeMatrixS_comps, method = "eucl", diag = F)

# And plot for visualization
aitHeatmap(aitchisonPlotB)



## ----- Vol Testing | Replicates ------------------------------
# Need to create a community data object--has R1 and R2 as separate samples.
any(rowSums(volumeMatrix) <1) # should be false
summary(colSums(volumeMatrix > 0)) # most ASVs only found on one site
summary(colSums(volumeMatrix))
hist(colSums(volumeMatrix))

#Create a curve within each letter sample for the three replicates. Basically
#trying to figure out how many new ATVs does each replicate contribute?

rarecurve(volumeMatrix[ grepl(pattern = "b-R1-A1|b-R1-A2|b-R1-A3",
                              x = rownames(volumeMatrix)), ])
rarecurve(volumeMatrix[ grepl(pattern = "b-R2-A1|b-R2-A2|b-R2-A3",
                                        x = rownames(volumeMatrix)), ])

# create species accumulation curves for each letter for buffers.
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-B1|b-R1-B2|b-R1-B3",
                                   x = rownames(volumeMatrix)), ]),
     xlab = "Number of replicates")
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-A1|b-R1-A2|b-R1-A3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 2)
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-C1|b-R1-C2|b-R1-C3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 3)
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-D1|b-R1-D2|b-R1-D3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 4)
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-E1|b-R1-E2|b-R1-E3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 5)
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-F1|b-R1-F2|b-R1-F3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 6)
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-G1|b-R1-G2|b-R1-G3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 7)
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-H1|b-R1-H2|b-R1-H3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 8)
plot(specaccum(volumeMatrix[ grepl(pattern = "b-R1-I1|b-R1-I2|b-R1-I3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 9)


# create species accumulation curves for both buffer and silica (pretending 6
# replicates instead of 3)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-B1|R1-B2|R1-B3",
                                   x = rownames(volumeMatrix)), ]),
     xlab = "Number of replicates [R1]", ylim = c(10,200))
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-A1|R1-A2|R1-A3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 2)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-C1|R1-C2|R1-C3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 3)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-D1|R1-D2|R1-D3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 4)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-E1|R1-E2|R1-E3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 5)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-F1|R1-F2|R1-F3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 6)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-G1|R1-G2|R1-G3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 7)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-H1|R1-H2|R1-H3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 8)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-I1|R1-I2|R1-I3",
                                   x = rownames(volumeMatrix)), ]),
     add = TRUE, col = 9)




# create species accumulation curves for both buffer and silica (pretending 6
# replicates instead of 3); remove rare species
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-B1|R1-B2|R1-B3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     xlab = "Number of replicates [R1, rare = 50]", ylim = c(10,120))
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-A1|R1-A2|R1-A3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 2)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-C1|R1-C2|R1-C3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 3)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-D1|R1-D2|R1-D3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 4)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-E1|R1-E2|R1-E3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 5)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-F1|R1-F2|R1-F3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 6)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-G1|R1-G2|R1-G3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 7)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-H1|R1-H2|R1-H3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 8)
plot(specaccum(volumeMatrix[ grepl(pattern = "R1-I1|R1-I2|R1-I3",
                                   x = rownames(volumeMatrix)), 
                             colSums(volumeMatrix) > rare]),
     add = TRUE, col = 9)





## ----- Vol Testing | Total Volume ----------------------------

plot(specaccum(volumeMatrixLetter[ grepl(pattern = "s-R1",
                                   x = rownames(volumeMatrixLetter)), ]),
     xlab = "Samples [silica & buffer, R1]")
plot(specaccum(volumeMatrixLetter[ grepl(pattern = "b-R1",
                                   x = rownames(volumeMatrixLetter)), ]),
     col = 2, add = T)




plot(specaccum(volumeMatrixLetter[ grepl(pattern = "R1",
                                   x = rownames(volumeMatrixLetter)), ]),
     xlab = "Samples [all R1]")

# what happens if we remove rare species?

plot(specaccum(volumeMatrixLetter[ grepl(pattern = "s-R1",
                                         x = rownames(volumeMatrixLetter)), 
                                   colSums(volumeMatrixLetter) > rare]),
     xlab = "Samples [silica & buffer, R1, rare = 50]")
plot(specaccum(volumeMatrixLetter[ grepl(pattern = "b-R1",
                                         x = rownames(volumeMatrixLetter)), 
                                   colSums(volumeMatrixLetter) > rare]),
     col = 2, add = T)
abline(h = 205)
abline(h = (205*.9), col = "blue")


plot(specaccum(volumeMatrixLetter[ grepl(pattern = "R1",
                                         x = rownames(volumeMatrixLetter)), 
                                   colSums(volumeMatrixLetter) > rare]),
     xlab = "Samples [all R1, rare = 50]")
abline(h = 205)
abline(h = (205*.9), col = "blue")
abline(h = (205*.75), col = "red")




rarecurve(volumeMatrixLetter[ grepl(pattern = "b-R1",
                              x = rownames(volumeMatrixLetter)), ])



iNEXTB <-
    iNEXT(speciesSite(volumeMatrixLetter[grepl(x = rownames(volumeMatrixLetter),
                                               pattern = "R1"),]),
          q = 0)
iNEXTB$DataInfo
iNEXTB$AsyEst
ggiNEXT(iNEXTB, type = 2)




## ----- CODE END ----------------------------------------------