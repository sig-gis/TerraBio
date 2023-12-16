# This analysis is based on soil samples collected at Reforesterra and processed by
# EcoMol. There are four land uses: CF = counter factual, Co = control/forest,
# Re = restauracao, Sy = syntropic. The 2023 analysis is the first year.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis.

# Code written Dec. 2023 by Karen Dyson

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
source("functions.R")

# script variable definitions
minlibrarySize = 5000
minRelativeAbund = 0.05
minAbsoluteAbund = 5

#Reitmeier, S., Hitch, T.C., Treichel, N., Fikas, N., Hausmann, B., Ramer-Tait,
#A.E., Neuhaus, K., Berry, D., Haller, D., Lagkouvardos, I. and Clavel, T.,
#2021. Handling of spurious sequences affects the outcome of high-throughput 16S
#rRNA gene amplicon profiling. ISME Communications, 1(1), pp.1-12.

phylum = c("Arthropoda") 
#phylum = c("Annelida", "Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca") #worms and insects



## Ingest 2023 data.

## common data
lookupColnames <- read.csv("lookupColnames.csv")
lookupSitenames2023 <- read.csv("Reforesterra/YR1_2023/lookupSitename2023.csv", strip.white = T)

infoColnames2023 <- c("N", "sampleID", "siteType", "labVolume",
                      "amplificationSuccess", "replicate",
                      "EcoMolID", "concentrationDNA_nguL",
                      "purityDNA")



# 2023 eDNA data import
        temp <- lookupColnames[order(lookupColnames$EM_OrderRTerra2023), ] 
        temp <- temp$TB_ColName[which(!is.na(temp$EM_ColNameRTerra2023) & !is.na(temp$EM_OrderRTerra2023))]

    rterra2023Raw <- read.csv("Reforesterra/YR1_2023/Reforesterra-Complete_analysis_results-2023-10-25.csv",
                             stringsAsFactors = F,
                             col.names = temp)
    
    
    ## WAITING ON THIS INFORMATION
    horta2023Info <- read.csv("Horta/Yr2_2023/HortaSamples2023.csv",
                              stringsAsFactors = F,
                              col.names = infoColnames2023,
                              header = T)
    

    ecoMolRawSummary(horta2022Raw)
    ecoMolRawSummary(horta2023Raw)
    
    horta2023Raw %>%
        group_by(identificationMaxTaxon) %>%
        summarise(asvAbsoluteAbundance)
    
    
    
## ----- Data cleaning and setup | 2022 -------------------------------

# Fix an error where EM91b-S1 should be EM91b-R2-Sy1
horta2022Raw$sample[horta2022Raw$sample == "EM91b-Sy1"] <- "EM91b-R2-Sy1"

# Remove the ASVHeader ">" character
horta2022Raw$ASVHeader <- str_sub(horta2022Raw$ASVHeader, 2, -1)

# Remove all of the volume testing information
horta2022Raw <- horta2022Raw[ horta2022Raw$metadata_4 %in% str_sub(horta2022Info$EcoMolID, -3, -1), ]

# Remove the R2 primer set information
horta2022Raw <- horta2022Raw[ horta2022Raw$primerName == "R1", ]

# Remove columns we don't need for analysis.
head(horta2022Raw)
horta2022Filt <- horta2022Raw[ , lookupColnames$TB_ColName[which(!is.na(lookupColnames$EM_ColName2022) &
                                                                    lookupColnames$Keep == "Y")]]
head(horta2022Filt)



# Look at phylum in the raw data
table(horta2022Filt$phylumBLASTn[horta2022Filt$preservation == "buffer"])

ggplot(horta2022Filt[horta2022Filt$preservation == "buffer",], aes(x = phylumBLASTn)) + 
    geom_bar()



# Create a group column
    horta2022Raw$LCGroupRep <- str_sub(horta2022Raw$metadata_4, 1, 3) %>% 
        str_replace_all(c("Co" = "Forest_",  "CF" = "Counterfactual_", "Re" = "Restoration_", "Sy" = "Syntropic_"))
    horta2022Raw$LCGroup <- str_split_fixed(horta2022Raw$LCGroupRep, "_",2)[,1] 
    
    horta2022Filt$LCGroupRep <- str_sub(horta2022Filt$metadata_4, 1, 3) %>% 
        str_replace_all(c("Co" = "Forest_",  "CF" = "Counterfactual_", "Re" = "Restoration_", "Sy" = "Syntropic_"))
    horta2022Filt$LCGroup <- str_split_fixed(horta2022Filt$LCGroupRep, "_",2)[,1] 


# Clean the data based on our quality variables
    # Remove sites not meeting minimum library size
    print(unique(horta2022Filt$sampleTotalAbd))
    removedSites <- unique(horta2022Filt$sample[horta2022Filt$sampleTotalAbd <= minlibrarySize])
    horta2022Filt <- horta2022Filt[ !(horta2022Filt$sample %in% removedSites) , ]
    stopifnot(length(unique(horta2022Filt$sample[horta2022Filt$sampleTotalAbd <= minlibrarySize])) == 0)

    print(removedSites)
    remove(removedSites)
    
    # Remove ASVs not meeting primer expected length
    horta2022Filt <- horta2022Filt[ horta2022Filt$primerExpectedLength == "in range", ]
    
    # Filter on the phylum
    horta2022Filt <- horta2022Filt[ horta2022Filt$phylumBLASTn %in% phylum, ]
    
    
    
# Diagnostic plot for buffer vs. silica 
    
    ggplot(data = horta2022Filt, 
           aes(x = preservation, y = sampleTotalAbd)) +
        geom_boxplot() + geom_point(aes(color = LCGroupRep))
    
    ggplot(data = horta2022Filt, 
           aes(asvAbsoluteAbundance)) +
        geom_density(aes(color = preservation)) +
        xlim(0,100)

    
# Either buffer or silica. Use silica preserved samples so 2022 and 2023 are the
# same. Use buffer samples to be consistent with 2022 analysis and because the
# results were a bit better...
    horta2022Filt <- horta2022Filt[ horta2022Filt$preservation == preservation2022, ]
    
    
    
## ----- Data cleaning and setup | 2023 -------------------------------
    
# Remove the ASVHeader ">" character
horta2023Raw$ASVHeader <- str_sub(horta2023Raw$ASVHeader, 2, -1)

# Remove columns we don't need for analysis.
head(horta2023Raw)
horta2023Filt <- horta2023Raw[ , lookupColnames$TB_ColName[which(lookupColnames$TB_ColName %in% temp &
                                                                    lookupColnames$Keep == "Y")]]
head(horta2023Filt)

# Create a grouping column
horta2023Raw$LCGroupRep <- str_sub(horta2023Raw$metadata_2, 1,2) %>% 
    str_replace_all(c("F" = "Forest_", "R" = "Restoration_", "S" = "Syntropic_", "V" = "Counterfactual_"))
horta2023Raw$LCGroup <- str_sub(horta2023Raw$metadata_2, 1,1) %>% 
    str_replace_all(c("F" = "Forest", "R" = "Restoration", "S" = "Syntropic", "V" = "Counterfactual"))

horta2023Filt$LCGroupRep <- str_sub(horta2023Filt$metadata_2, 1,2) %>% 
    str_replace_all(c("F" = "Forest_", "R" = "Restoration_", "S" = "Syntropic_", "V" = "Counterfactual_"))
horta2023Filt$LCGroup <- str_sub(horta2023Filt$metadata_2, 1,1) %>% 
    str_replace_all(c("F" = "Forest", "R" = "Restoration", "S" = "Syntropic", "V" = "Counterfactual"))

table(horta2023Filt$phylumBLASTn)

    # Clean the data based on our quality variables
    # Remove sites not meeting minimum library size -- none fail
    removedSites <- unique(horta2023Filt$sample[horta2023Filt$sampleTotalAbd <= minlibrarySize])
    horta2023Filt <- horta2023Filt[ !(horta2023Filt$sample %in% removedSites) , ]
    stopifnot(length(unique(horta2023Filt$sample[horta2023Filt$sampleTotalAbd <= minlibrarySize])) == 0)
    
    print(removedSites)
    remove(removedSites)
    
    # Remove ASVs not meeting primer expected length
    horta2023Filt <- horta2023Filt[ horta2023Filt$primerExpectedLength == "in range", ]
    
    # Filter on the phylum
    horta2023Filt <- horta2023Filt[ horta2023Filt$phylumBLASTn %in% phylum, ]
    
        remove(temp)
        
        
        ecoMolFiltSummary(horta2022Filt)
        ecoMolFiltSummary(horta2023Filt) # many more reads after filtering
        
        
# This identifies sequences in both years' data--this suggests it is probably a
# "real" sequence and not a fluke. It also facilitates the joint PCAs
    commonASV <- tibble(
        ASVSequence = horta2023Raw$ASVSequence[horta2023Raw$ASVSequence %in% horta2022Raw$ASVSequence]
    )
    
    commonASV <- merge(x = commonASV,
                       y = horta2022Raw[ , colnames(horta2022Raw) %in% c("ASVSequence", "ASVHeader")],
                       by = "ASVSequence",
                       all.x = T
                       ) %>% unique()
    colnames(commonASV)[2] <- "ASVHeader2022"
    commonASV <- merge(x = commonASV,
                       y = horta2023Raw[ , colnames(horta2023Raw) %in% c("ASVSequence", "ASVHeader")],
                       by = "ASVSequence",
                       all.x = T
    ) %>% unique()
    colnames(commonASV)[3] <- "ASVHeader2023"
    commonASV$lookup <- paste0(commonASV$ASVHeader2022, "=",commonASV$ASVHeader2023)


# ----- Data quality checks | 2022 & 2023 -------------------------------

# Compare 2022 and 2023 DNA concentration.

p0 <- ggplot(subset(horta2022Info, storage == "Silica"),
             aes(x = volumeSampleID, y = concentrationDNA_nguL)) + 
    geom_boxplot() + geom_point(color = "green") +
    ylim(5,80) + xlab("LC type (silica 2022)")
    
p1 <- ggplot(subset(horta2022Info, storage == "Buffer"),
             aes(x = volumeSampleID, y = concentrationDNA_nguL)) + 
             geom_boxplot() + geom_point(color = "darkgreen") +
        ylim(5,80) + xlab("LC type (buffer 2022)")
    
p2 <- ggplot(horta2023Info,
             aes(x = siteType, y = concentrationDNA_nguL)) + 
             geom_boxplot() + geom_point(aes(color = labVolume)) +
    ylim(5,80) + xlab("LC type (silica 2023)")
    

    grid.arrange(p0, p1, p2, nrow=2)
    
# Compare 2022 and 2023 DNA purity.
    
    p1 <- ggplot(subset(horta2022Info, storage == "Buffer"),
                 aes(x = volumeSampleID, y = purityDNA)) + 
        geom_boxplot() + geom_point(color = "darkgreen") +
        ylim(1.45,2) + xlab("LC type (buffer 2022)")
    
    p2 <- ggplot(horta2023Info,
                 aes(x = siteType, y = purityDNA)) + 
        geom_boxplot() + geom_point(aes(color = labVolume)) +
        ylim(1.45,2) + xlab("LC type (silica 2023)")
    
    
    grid.arrange(p1, p2, ncol=2)    
    
    remove(p0, p1, p2)
    
# overall silica data quality of 2023 appears better than 2022, with the 2023
# data having higher concentration and purity. For the concentration, forest and
# restoration had low concentration in 2022 while forest and counterfactual had
# the lowest in 2023. the forest had low volume sent to the lab and the first
# amplification failed. Purity for the forest was fine, right in the middle of
# the pack. however we might expect the forest to have concentrations more
# similar to the Restoration sites. final note, for 2022 the buffer had much
# higher concentration than the silica samples.
    
# Look at ASV absolute abundance
    
    ggplot(horta2022Raw, aes(x = LCGroup, y = asvAbsoluteAbundance)) + 
        geom_boxplot() 
    ggplot(horta2022Filt, aes(x = LCGroup, y = asvAbsoluteAbundance)) + 
        geom_boxplot() 
    
    ggplot(horta2023Raw, aes(x = LCGroup, y = asvAbsoluteAbundance)) + 
        geom_boxplot() #+ ylim(0,10000)
    ggplot(horta2023Filt, aes(x = LCGroup, y = asvAbsoluteAbundance)) + 
        geom_boxplot() #+ ylim(0,10000)
    # A few very high abundance ASVs in Forest and Restoration. many fewer for
    # Syntropic and the Pasture, except for one Pasture site.
    
    #what are the weird ASVs?
    horta2023Raw[horta2023Raw$asvAbsoluteAbundance > 5000, ]
    horta2023Filt[horta2023Filt$asvAbsoluteAbundance > 5000, ]
    # One Coleoptera in the Restoration, one Anoplotermes (termite) in the
    # Forest. The raw data also has an Annelid Metaphire sp.
    
# look at sample total abundance
    ggplot(horta2022Filt, aes(x = LCGroup, y = sampleTotalAbd)) + 
        geom_point(aes(color = LCGroupRep))
    
    ggplot(horta2023Filt, aes(x = LCGroup, y = sampleTotalAbd)) + 
        geom_point(aes(color = LCGroupRep)) 
    
# what is the one weird V?
    horta2023Raw[horta2023Raw$sampleTotalAbd > 50000, ]
    sum(horta2023Raw$asvAbsoluteAbundance[horta2023Raw$metadata_2 == "V2"])
    # things that get filtered out in the process... none of which are identified.
    horta2023Raw[ horta2023Raw$asvAbsoluteAbundance > 750 & horta2023Raw$metadata_2 == "V2", ]
 
    
    
    
# ----- Filter for rare species | 2022 & 2023 -------------------------

# We want to remove rare species that could be errors 

# Remove ASVs below minimum relative abundance on sample. 
    plot(
        horta2022Filt$asvAbsoluteAbundance,
        horta2022Filt$relativeAbundanceOnSample,
        ylim = c(0, .05),
        xlim = c(0, 50)
    )
    abline(v = minAbsoluteAbund, col = "red")
    
    print("2022 filt. ASV < min absolute abundance:")
    print(table(horta2022Filt$asvAbsoluteAbundance < minAbsoluteAbund))
    
    horta2022Filt <- horta2022Filt[ horta2022Filt$asvAbsoluteAbundance >= minAbsoluteAbund, ]
    
# 2023: Remove ASVs below minimum relative abundance on sample. 
    plot(
        horta2023Filt$asvAbsoluteAbundance,
        horta2023Filt$relativeAbundanceOnSample,
        ylim = c(0, .003),
        xlim = c(0, 50)
    )
    abline(v = minAbsoluteAbund, col = "red")
    
    # 2023 has very low relative abundance on sample compared with Horta 2022 and Cafe Apui 2022
    
    print("2023 filt. ASV < min absolute abundance:")
    print(table(horta2023Filt$asvAbsoluteAbundance < minAbsoluteAbund))

    horta2023Filt <- horta2023Filt[ horta2023Filt$asvAbsoluteAbundance >= minAbsoluteAbund, ]
    
    
# ----- Create Matrices | 2022 ---------------------------------------
    
    # Create a matrix with replicates as individual "sites"
    
    hortaMatrix2022 <- ez.matrify(horta2022Filt, species.name = "ASVHeader",
                                  site.name = "sample", abundance = "asvAbsoluteAbundance")
    
    hist(colSums(hortaMatrix2022), breaks = 50)
    # Create a matrix where the replicates for land use types are combined
    
    hortaType2022 <- horta2022Filt %>%
        dplyr::select(LCGroup, ASVHeader, asvAbsoluteAbundance) %>%
        group_by(LCGroup, ASVHeader) %>%
        summarise(abundance = sum(asvAbsoluteAbundance))
    
    hortaMatrixType2022 <- ez.matrify(hortaType2022, species.name = "ASVHeader",
                                      site.name = "LCGroup", abundance = "abundance")
    
    #test to make sure everything got in
    any((colSums(hortaMatrixType2022)-colSums(hortaMatrix2022)) != 0 )
    
# ----- Create Matrices | 2023 ---------------------------------------
    
    # Create a matrix with replicates as individual "sites"
    
    hortaMatrix2023 <- ez.matrify(horta2023Filt, species.name = "ASVHeader",
                                  site.name = "sample", abundance = "asvAbsoluteAbundance")
    
    hist(colSums(hortaMatrix2023), breaks = 100)
    
    
    
    
    
    # Create a matrix where the replicates for land use types are combined
    
    hortaType2023 <- horta2023Filt %>%
        dplyr::select(LCGroup, ASVHeader, asvAbsoluteAbundance) %>%
        group_by(LCGroup, ASVHeader) %>%
        summarise(abundance = sum(asvAbsoluteAbundance))
    
    hortaMatrixType2023 <- ez.matrify(hortaType2023, species.name = "ASVHeader",
                                      site.name = "LCGroup", abundance = "abundance")
    
    #test to make sure everything got in
    any((colSums(hortaMatrixType2023)-colSums(hortaMatrix2023)) != 0 )
    
    
    