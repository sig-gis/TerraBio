# This analysis is based on soil samples collected at Horta and processed by
# EcoMol. There are four land uses: CF = counter factual, Co = control/forest,
# Re = restauracao, Sy = syntropic. The 2023 analysis combines both 2022 and
# 2023 data.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis.

# Code written Sept. 2023 by Karen Dyson

# to do list:
#   + need to look at which ASVs should be combined or not...
#   + check high ASV loss. approx 1000-1200 ASV rows are kept out of around 20k.
#   + add diagnostics/plotting
#   + add diagnostics for what ASVs are lost


## ----- Data ingestion ---------------------------------------

# libraries
library(ggplot2)
library(gridExtra)
library(tidyr)
library(vegan)
library(stringr)
library(dplyr)

# Ingest codes
source("../../../RCode/R_Scripts/triplet_fixer.R")

# script variable definitions
minlibrarySize = 4000
minRelativeAbund = 0.05
minAbsoluteAbund = 5

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
    


# ----- Data quality checks | 2022 & 2023 -------------------------------

# Compare 2022 and 2023 DNA concentration.

p1 <- ggplot(subset(horta2022Info, storage == "Silica"),
             aes(x = volumeSampleID, y = concentrationDNA_nguL)) + 
             geom_boxplot() + geom_point(color = "darkgreen") +
        ylim(5,40)
    
p2 <- ggplot(horta2023Info,
             aes(x = siteType, y = concentrationDNA_nguL)) + 
             geom_boxplot() + geom_point(aes(color = labVolume)) +
    ylim(5,40)
    

    grid.arrange(p1, p2, ncol=2)
    
# Compare 2022 and 2023 DNA purity.
    
    p1 <- ggplot(subset(horta2022Info, storage == "Silica"),
                 aes(x = volumeSampleID, y = purityDNA)) + 
        geom_boxplot() + geom_point(color = "darkgreen") +
        ylim(1.45,2)
    
    p2 <- ggplot(horta2023Info,
                 aes(x = siteType, y = purityDNA)) + 
        geom_boxplot() + geom_point(aes(color = labVolume)) +
        ylim(1.45,2)
    
    
    grid.arrange(p1, p2, ncol=2)    
    
# overall data quality of 2023 appears better than 2022, with the 2023 data
# having higher concentration and purity. For the concentration, forest and
# restoration had low concentration in 2022 while forest and counterfactual had
# the lowest in 2023. the forest had low volume sent to the lab and the first
# amplification failed. Purity for these was fine also, right in the middle of
# the pack. final note, for 2022 the buffer had much higher concentration than
# the silica samples.
    
# Look at ASV absolute abundance
    
    ggplot(horta2022Raw, aes(x = metadata_4, y = asvAbsoluteAbundance)) + 
        geom_boxplot() #+ geom_jitter()
    
    ggplot(horta2023Raw, aes(x = metadata_4, y = asvAbsoluteAbundance)) + 
        geom_boxplot() #+ geom_jitter()
    
    
    # look at sample total abundance
    #
    
    
    ## ----- Storage method graphs --------------------------------
    
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
    
    