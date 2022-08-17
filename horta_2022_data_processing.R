# This analyis is based on soil samples collected in June at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis (Both volume testing and Horta TB indicators)

## ----- Data ingestion ---------------------------------------

lookupColnames <- read.csv("Horta/lookupColnames.csv")
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
                      col.names = infoColnames)

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

rare = 50


## ----- Data cleaning and setup -------------------------------
# Remove columns we don't need for analysis.
head(hortaAllASV)
hortaAllASV <- hortaAllASV[ , lookupColnames$Keep == "Y"]
head(hortaAllASV)

# Remove ASVs that don't meet criteria.

hortaAllASV <- hortaAllASV[ hortaAllASV$primerExpectedLength == "in range", ]

table(hortaAllASV$phylumBLASTn[grepl(pattern = "R2", x = hortaAllASV$sample)])
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

hortaAllASV <- hortaAllASV[ hortaAllASV$phylumBLASTn == "Arthropoda", ]

# Now we need to address sample names and split out the volume samples.
table(hortaAllASV$sample)
# A1, A2, A3 etc. A is the volume sample (9 total) and 1-3 is the replicate.
# CF = counter factual, Co = control/forest, Re = restauracao, Sy = syntropic.
# b = buffer, s = silica.
# R1 and R2 are the two primer sets.

hortaString <- c("CF", "Co", "Re", "Sy")

# Split the horta and volume samples.
hortaASV <- hortaAllASV[ grepl(pattern = paste0(hortaString, collapse = "|") , x = hortaAllASV$metadata_4) , ]
volumeASV <- hortaAllASV[ !(hortaAllASV$metadata_4 %in% hortaASV$metadata_4) , ]
volumeASV$ASVHeader <- str_sub(volumeASV$ASVHeader, 2, -1)

volumeMatrix <- ez.matrify(volumeASV, species.name = "ASVHeader",
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


remove(hortaAllASV, hortaString)
