# This analysis is based on soil samples collected at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis.

# Code written Sept. 2023 by Karen Dyson


## ----- Data ingestion ---------------------------------------

# libraries
library(ggplot2)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(dplyr)

# Ingest codes
source("../allianceBranding.R")
source("../functions.R")
source("../multiyear_functions.R")
source("../../../RCode/R_Scripts/triplet_fixer.R")

# script variable definitions
minlibrarySize = 5000
minRelativeAbund = 0.05
minAbsoluteAbund = 5
rare = 50

#Reitmeier, S., Hitch, T.C., Treichel, N., Fikas, N., Hausmann, B., Ramer-Tait,
#A.E., Neuhaus, K., Berry, D., Haller, D., Lagkouvardos, I. and Clavel, T.,
#2021. Handling of spurious sequences affects the outcome of high-throughput 16S
#rRNA gene amplicon profiling. ISME Communications, 1(1), pp.1-12.

phylum = c("Arthropoda")
#phylum = c("Annelida", "Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca") #worms and insects



## First ingest 2022 and 2023 data.

## common data
lookupColnames <- read.csv("lookupColnames.csv")
lookupSitenames <- read.csv("lookupSitename.csv")

infoColnames <- c("N", "storage", "volumeSampleID",
                  "replicate", "replicateVolume",
                  "EcoMolID", "concentrationDNA_nguL",
                  "purityDNA")

# 2022 data--unique case because we had the volume information to sort out
horta2022Raw <- read.csv("Horta_2022_bioinfo_results.csv",
                        stringsAsFactors = F,
                        col.names = lookupColnames$TB_ColName)
volumeInfo <- read.csv("VolumeTest2022.csv",
                       stringsAsFactors = F,
                       col.names = infoColnames)
horta2022Info <- read.csv("HortaSamples2022.csv",
                      stringsAsFactors = F,
                      col.names = infoColnames, header = T)

# 2023 data
horta2023Raw <- read.csv("Horta_2023_bioinfo_results.csv",
                         stringsAsFactors = F,
                         col.names = lookupColnames$TB_ColName)
horta2023Info <- read.csv("HortaSamples2023.csv",
                          stringsAsFactors = F,
                          col.names = infoColnames, header = T)


## ----- Data cleaning and setup | 2022 -------------------------------

# Fix an error where EM91b-S1 should be EM91b-R2-Sy1
hortaAllASV$sample[hortaAllASV$sample == "EM91b-Sy1"] <- "EM91b-R2-Sy1"


# Remove columns we don't need for analysis.
head(horta2022Raw)
horta2022Raw <- horta2022Raw[ , lookupColnames$Keep == "Y"]
head(horta2022Raw)



# to do list:
#   + add column names for 2023
#   + import 2023 data (even if wrong) to build out
#   + need to look at which ASVs should be combined or not...






