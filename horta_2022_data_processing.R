# This analysis is based on soil samples collected in June at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis (Both volume testing and Horta TB indicators)

## ----- Data ingestion ---------------------------------------

lookupColnames <- read.csv("Horta/lookupColnames.csv")
lookupSitenames <- read.csv("Horta/lookupSitename.csv")

infoColnames <- c("N", "storage", "volumeSampleID",
                  "replicate", "replicateVolume",
                  "EcoMolID", "concentrationDNA_nguL",
                  "purityDNA")

hortaAllASV <- read.csv("Horta/Horta_2022_bioinfo_results.csv",
                        stringsAsFactors = F,
                        col.names = lookupColnames$TB_ColName)
volumeInfo <- read.csv("Horta/VolumeTest.csv",
                       stringsAsFactors = F,
                       col.names = infoColnames)
hortaInfo <- read.csv("Horta/HortaSamples.csv",
                      stringsAsFactors = F,
                      col.names = infoColnames, header = T)



# Fix an error where EM91b-S1 should be EM91b-R2-Sy1
hortaAllASV$sample[hortaAllASV$sample == "EM91b-Sy1"] <- "EM91b-R2-Sy1"

# libraries
library(ggplot2)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(dplyr)


source("allianceBranding.R")

source("functions.R")
source("../../../RCode/R_Scripts/triplet_fixer.R")

# script variable definitions
minlibrarySize = 5000
minRelativeAbund = 0.05
minAbsoluteAbund = 5

#phylum = c("Annelida", "Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca") #worms and insects
#phylum = c("Rhodophyta", "Streptophyta", "Chlorophyta", "Bacillariophyta") #plants and algae
#phylum = c("Rotifera", "Tubulinea", "Oomycota", "Gastrotricha", "Discosea", "Endomyxa", "Ascomycota") #microorganisms and fungi
phylum = c("Arthropoda")




## ----- Data cleaning and setup -------------------------------
# Remove columns we don't need for analysis.
head(hortaAllASV)
hortaAllASV <- hortaAllASV[ , lookupColnames$Keep == "Y"]
head(hortaAllASV)

# Remove ASVs that don't meet criteria.

hortaAllASV <- hortaAllASV[ hortaAllASV$primerExpectedLength == "in range", ]

    # EM91b-R2-Re2 disappears after filtering, find out why.
    table(hortaAllASV$phylumBLASTn[grepl(pattern = "b-R2-Re2", x = hortaAllASV$sample)])
    # answer--no Arthropods detected.

    # check sums
    
    table(hortaAllASV$phylumBLASTn[grepl(pattern = "R1", x = hortaAllASV$sample)])
    nrow(hortaAllASV[grepl(
        pattern = "R1", x = hortaAllASV$sample) & hortaAllASV$phylumBLASTn == "Arthropoda", ])
    # Both R1 and R2
    sum(hortaAllASV$phylumBLASTn == "Arthropoda")/
        (sum(table(hortaAllASV$phylumBLASTn))-44573)
    #R1 only
    sum(hortaAllASV$phylumBLASTn[grepl(
        pattern = "R1", x = hortaAllASV$sample)] == "Arthropoda")/
        (sum(table(hortaAllASV$phylumBLASTn[grepl(
            pattern = "R1", x = hortaAllASV$sample)]))-43224)
    
    #R2 only
    sum(hortaAllASV$phylumBLASTn[grepl(
        pattern = "R2", x = hortaAllASV$sample)] == "Arthropoda")/
        (sum(table(hortaAllASV$phylumBLASTn[grepl(
            pattern = "R2", x = hortaAllASV$sample)]))-1329)

hortaAllASV <- hortaAllASV[ hortaAllASV$phylumBLASTn %in% phylum, ]

# Now we need to address sample names and split out the volume samples.
table(hortaAllASV$sample)
# A1, A2, A3 etc. A is the volume sample (9 total) and 1-3 is the replicate.
# CF = counter factual, Co = control/forest, Re = restauracao, Sy = syntropic.
# b = buffer, s = silica.
# R1 and R2 are the two primer sets.

hortaString <- c("CF", "Co", "Re", "Sy")

# Split the horta and volume samples.
hortaAllASV$ASVHeader <- str_sub(hortaAllASV$ASVHeader, 2, -1)

hortaASV <- hortaAllASV[ grepl(pattern = paste0(hortaString, collapse = "|") , x = hortaAllASV$metadata_4) , ]
volumeASV <- hortaAllASV[ !(hortaAllASV$metadata_4 %in% hortaASV$metadata_4) , ]


volumeMatrix <- ez.matrify(volumeASV, species.name = "ASVHeader",
                           site.name = "sample", abundance = "asvAbsoluteAbundance")

hortaMatrix <- ez.matrify(hortaASV, species.name = "ASVHeader",
                          site.name = "sample", abundance = "asvAbsoluteAbundance")

# For the volume samples, create a dataset where the samples are combined by
# replicate.

volumeASV$sampleLetter <- paste0(volumeASV$metadata_5, "-", volumeASV$primerName, "-", str_sub(volumeASV$metadata_4, 1, 1))

abundanceLetter <- volumeASV %>%
    dplyr::select(sampleLetter, ASVHeader, asvAbsoluteAbundance) %>%
    group_by(sampleLetter, ASVHeader) %>%
    summarise(abundance = sum(asvAbsoluteAbundance))

volumeMatrixLetter <- ez.matrify(abundanceLetter, species.name = "ASVHeader",
                                 site.name = "sampleLetter", abundance = "abundance")

#test to make sure everything got in
any((colSums(volumeMatrixLetter)-colSums(volumeMatrix)) > 0 )

# For the horta samples, create a dataset where the land use type are combined
# by replicate.

hortaASV$sampleLetter <- paste0(hortaASV$metadata_5, "-", hortaASV$primerName, "-", str_sub(hortaASV$metadata_4, 1, 2))

hortaLetter <- hortaASV %>%
    dplyr::select(sampleLetter, ASVHeader, asvAbsoluteAbundance) %>%
    group_by(sampleLetter, ASVHeader) %>%
    summarise(abundance = sum(asvAbsoluteAbundance))

hortaMatrixLetter <- ez.matrify(hortaLetter, species.name = "ASVHeader",
                                 site.name = "sampleLetter", abundance = "abundance")

#test to make sure everything got in
any((colSums(hortaMatrixLetter)-colSums(hortaMatrix)) > 0 )

remove(hortaAllASV)
