# This analysis is based on soil samples collected at Horta and processed
# by EcoMol. The 2023 analysis combines both 2022 and 2023 data.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis.

# Code written Sept. 2023 by Karen Dyson

# to do list:
#   + need to look at which ASVs should be combined or not...
#   + check high ASV loss. approx 1000-1200 ASV rows are kept out of around 20k.
#   + add diagnostics/plotting


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
source("allianceBranding.R")
source("functions.R")
source("multiyear_functions.R")
source("../../../RCode/R_Scripts/triplet_fixer.R")

# script variable definitions
minlibrarySize = 4000
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
lookupSitenames2022 <- read.csv("Horta/YR1_2022/lookupSitename2022.csv")
lookupSitenames2023 <- read.csv("Horta/YR2_2023/lookupSitename2023.csv")

infoColnames2022 <- c("N", "storage", "volumeSampleID",
                  "replicate", "replicateVolume",
                  "EcoMolID", "concentrationDNA_nguL",
                  "purityDNA")
infoColnames2023 <- c("N", "sampleID", "siteType", "labVolume",
                      "amplificationSuccess", "replicate",
                      "EcoMolID", "concentrationDNA_nguL",
                      "purityDNA")

# 2022 eDNA data import

# CF = counter factual, Co = control/forest, Re = restauracao, Sy = syntropic.

# Note that 2022 was a unique case because we had the volume information, primer
# information, and preservation method information to sort out
        temp <- lookupColnames$TB_ColName[which(!is.na(lookupColnames$EM_ColName2022))]
        
    horta2022Raw <- read.csv("Horta/Yr1_2022/Horta_2022_bioinfo_results.csv",
                            stringsAsFactors = F,
                            col.names = temp) 
    horta2022Info <- read.csv("Horta/Yr1_2022/HortaSamples2022.csv",
                          stringsAsFactors = F,
                          col.names = infoColnames2022,
                          header = T)

# 2023 eDNA data import
        temp <- lookupColnames[order(lookupColnames$EM_Order2023), ] 
        temp <- temp$TB_ColName[which(!is.na(temp$EM_ColName2023) & !is.na(temp$EM_Order2023))]

    horta2023Raw <- read.csv("Horta/Yr2_2023/Horta_2023_bioinfo_results.csv",
                             stringsAsFactors = F,
                             col.names = temp)
    horta2023Info <- read.csv("Horta/Yr2_2023/HortaSamples2023.csv",
                              stringsAsFactors = F,
                              col.names = infoColnames2023,
                              header = T)
    


## ----- Data cleaning and setup | 2022 -------------------------------

# Fix an error where EM91b-S1 should be EM91b-R2-Sy1
horta2022Raw$sample[horta2022Raw$sample == "EM91b-Sy1"] <- "EM91b-R2-Sy1"

# Remove the ASVHeader ">" character
horta2022Raw$ASVHeader <- str_sub(horta2022Raw$ASVHeader, 2, -1)

# Remove columns we don't need for analysis.
head(horta2022Raw)
horta2022Raw <- horta2022Raw[ , lookupColnames$TB_ColName[which(!is.na(lookupColnames$EM_ColName2022) & lookupColnames$Keep == "Y")]]
head(horta2022Raw)

# Remove all of the volume testing information
horta2022Raw <- horta2022Raw[ horta2022Raw$metadata_4 %in% str_sub(horta2022Info$EcoMolID, -3, -1), ]

# Remove the R2 primer set information
horta2022Raw <- horta2022Raw[ horta2022Raw$primerName == "R1", ]

# Use just the silica preserved samples so 2022 and 2023 are the same. I think
# we used the buffer for the 2022 data because the results were a bit better...
horta2022Raw <- horta2022Raw[ horta2022Raw$preservation == "silica", ]

# Clean the data based on our quality variables
    # Remove sites not meeting minimum library size
    removedSites <- unique(horta2022Raw$sample[horta2022Raw$sampleTotalAbd <= minlibrarySize])
    horta2022Raw <- horta2022Raw[ !(horta2022Raw$sample %in% removedSites) , ]
    stopifnot(length(unique(horta2022Raw$sample[horta2022Raw$sampleTotalAbd <= minlibrarySize])) == 0)

    print(removedSites)
    remove(removedSites)
    
    # Remove ASVs not meeting primer expected length
    horta2022Raw <- horta2022Raw[ horta2022Raw$primerExpectedLength == "in range", ]
    
    # Filter on the phylum
    horta2022Raw <- horta2022Raw[ horta2022Raw$phylumBLASTn %in% phylum, ]
    

    
    
## ----- Data cleaning and setup | 2023 -------------------------------
    
# Remove the ASVHeader ">" character
horta2023Raw$ASVHeader <- str_sub(horta2023Raw$ASVHeader, 2, -1)

# Remove columns we don't need for analysis.
head(horta2023Raw)
horta2023Raw <- horta2023Raw[ , lookupColnames$TB_ColName[which(lookupColnames$TB_ColName %in% temp & lookupColnames$Keep == "Y")]]
head(horta2023Raw)

    # Clean the data based on our quality variables
    # Remove sites not meeting minimum library size -- none fail
    removedSites <- unique(horta2023Raw$sample[horta2023Raw$sampleTotalAbd <= minlibrarySize])
    horta2023Raw <- horta2023Raw[ !(horta2023Raw$sample %in% removedSites) , ]
    stopifnot(length(unique(horta2023Raw$sample[horta2023Raw$sampleTotalAbd <= minlibrarySize])) == 0)
    
    print(removedSites)
    remove(removedSites)
    
    # Remove ASVs not meeting primer expected length
    horta2023Raw <- horta2023Raw[ horta2023Raw$primerExpectedLength == "in range", ]
    
    # Filter on the phylum
    horta2023Raw <- horta2023Raw[ horta2023Raw$phylumBLASTn %in% phylum, ]
    
        remove(temp)
    
# ----- Create Matrices | 2022 ---------------------------------------
    
# Create a matrix with replicates as individual "sites"
    
    hortaMatrix2022 <- ez.matrify(horta2022Raw, species.name = "ASVHeader",
                              site.name = "sample", abundance = "asvAbsoluteAbundance")
    
        hist(colSums(hortaMatrix2022), breaks = 50)
# Create a matrix where the replicates for land use types are combined
    
    horta2022Raw$sampleLetter <- paste0(horta2022Raw$metadata_5, "-", horta2022Raw$primerName, "-", str_sub(horta2022Raw$metadata_4, 1, 2))
    
    hortaType2022 <- horta2022Raw %>%
        dplyr::select(sampleLetter, ASVHeader, asvAbsoluteAbundance) %>%
        group_by(sampleLetter, ASVHeader) %>%
        summarise(abundance = sum(asvAbsoluteAbundance))
    
    hortaMatrixType2022 <- ez.matrify(hortaType2022, species.name = "ASVHeader",
                                    site.name = "sampleLetter", abundance = "abundance")
    
#test to make sure everything got in
    any((colSums(hortaMatrixType2022)-colSums(hortaMatrix2022)) != 0 )
    
# ----- Create Matrices | 2022 ---------------------------------------
    
# Create a matrix with replicates as individual "sites"
    
    hortaMatrix2023 <- ez.matrify(horta2023Raw, species.name = "ASVHeader",
                                  site.name = "sample", abundance = "asvAbsoluteAbundance")
    
    hist(colSums(hortaMatrix2023), breaks = 100)
# Create a matrix where the replicates for land use types are combined
    
    horta2023Raw$sampleLetter <- str_sub(horta2023Raw$metadata_2, 1,1)

    hortaType2023 <- horta2023Raw %>%
        dplyr::select(sampleLetter, ASVHeader, asvAbsoluteAbundance) %>%
        group_by(sampleLetter, ASVHeader) %>%
        summarise(abundance = sum(asvAbsoluteAbundance))
    
    hortaMatrixType2023 <- ez.matrify(hortaType2023, species.name = "ASVHeader",
                                    site.name = "sampleLetter", abundance = "abundance")
    
#test to make sure everything got in
    any((colSums(hortaMatrixType2023)-colSums(hortaMatrix2023)) != 0 )
    



