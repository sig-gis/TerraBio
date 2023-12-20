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
    
    
    
## ----- Data cleaning and setup | 2023 -------------------------------

# Remove the ASVHeader ">" character
    rterra2023Raw$ASVHeader <- str_sub(rterra2023Raw$ASVHeader, 2, -1)

# Remove columns we don't need for analysis.
head(rterra2023Raw)
rterra2023Filt <- rterra2023Raw[ , lookupColnames$TB_ColName[which(!is.na(lookupColnames$EM_ColNameRTerra2023) &
                                                                    lookupColnames$Keep == "Y")]]
head(rterra2023Filt)

# Remove the controls in the filtered data
rterra2023Filt <- rterra2023Filt[!(rterra2023Filt$sample == "EM135c3-Neg"), ]

# Look at phylum in the raw data
table(rterra2023Filt$phylumBLASTn)

ggplot(rterra2023Filt, aes(x = phylumBLASTn)) + 
    geom_bar() + coord_flip()



# Create a group columns (first should be lc-site-replicates, second is lc-site, third is prettified lc only)

# For Reforesterra in 2023, metadata_2 is the site number, _3 through _5 are different codes for the lc, and _7 is the replicate. Can do this either by using the lookup table (merge + paste) or str_sub/str_replace_all. The challenge is that they have intervention A and B for one of the sites. There are also two types of counterfactual.



rterra2023Filt <- merge(rterra2023Filt,
                        lookupSitenames2023[ , 2:10],
                        by.x = "sample", by.y = "sample.code")

# lc-site-replicates and lc-site; lc only handled by "treatment"
rterra2023Filt <- rterra2023Filt %>%
    mutate(LCSiteRep = paste0(treatment.code, "-S", site.code, "-R", replicate),
           LCSite    = paste0(treatment.code, "-S", site.code))


# Clean the data based on our quality variables
    # Remove sites not meeting minimum library size
    print(unique(rterra2023Filt$sampleTotalAbd))
    rterra2023Filt %>% select(sample, sampleTotalAbd) %>% unique()
    removedSites <- unique(rterra2023Filt$sample[rterra2023Filt$sampleTotalAbd <= minlibrarySize])
    rterra2023Filt <- rterra2023Filt[ !(rterra2023Filt$sample %in% removedSites) , ]
    stopifnot(length(unique(rterra2023Filt$sample[rterra2023Filt$sampleTotalAbd <= minlibrarySize])) == 0)

    print(removedSites)
    remove(removedSites)
    
    # Remove ASVs not meeting primer expected length
    rterra2023Filt <- rterra2023Filt[ rterra2023Filt$primerExpectedLength == "in range", ]
    
    # Filter on the phylum
    rterra2023Filt <- rterra2023Filt[ rterra2023Filt$phylumBLASTn %in% phylum, ]
    
  

# ----- * Data quality checks | 2023 -------------------------------

# Examine 2023 DNA concentration.
    
    # Waiting on data from EcoMOL

p0 <- ggplot(subset(horta2022Info, storage == "Silica"),
             aes(x = volumeSampleID, y = concentrationDNA_nguL)) + 
    geom_boxplot() + geom_point(color = "green") +
    ylim(5,80) + xlab("LC type (silica 2022)")
    
    grid.arrange(p0, p1, p2, nrow=2)
    
# Examine 2023 DNA purity.
    
    p1 <- ggplot(subset(horta2022Info, storage == "Buffer"),
                 aes(x = volumeSampleID, y = purityDNA)) + 
        geom_boxplot() + geom_point(color = "darkgreen") +
        ylim(1.45,2) + xlab("LC type (buffer 2022)")
    

    grid.arrange(p0, p1, ncol=2)    
    
    remove(p0, p1, p2)
    
# add analysis notes here
    
# Look at ASV absolute abundance

    ggplot(rterra2023Filt, aes(x = treatment, y = asvAbsoluteAbundance)) + 
        geom_boxplot() 
    # if needed, repeat for Raw data.

    #what are the v. high ASVs?
    rterra2023Filt[rterra2023Filt$asvAbsoluteAbundance > 5000, ]
    # One Coleoptera (possibly lesser mealworm Alphitobius diaperinus or relative) found in Riparian pasture. 
    
# look at sample total abundance

    ggplot(rterra2023Filt, aes(x = LCSite, y = sampleTotalAbd)) + 
        geom_point() 
    

    
# ----- Filter for rare species | 2023 -------------------------

# We want to remove rare species that could be errors 

# Remove ASVs below minimum relative abundance on sample. 
    plot(
        rterra2023Filt$asvAbsoluteAbundance,
        rterra2023Filt$relativeAbundanceOnSample,
        ylim = c(0, .01),
        xlim = c(0, 50)
    )
    abline(v = minAbsoluteAbund, col = "red")
    
    print("2023 filt. ASV < min absolute abundance:")
    print(table(rterra2023Filt$asvAbsoluteAbundance < minAbsoluteAbund))
    
    rterra2023Filt <- rterra2023Filt[ rterra2023Filt$asvAbsoluteAbundance >= minAbsoluteAbund, ]
    

# ----- Create Matrices | 2023 ---------------------------------------
    
    # Create a matrix with replicates as individual "sites"
    
    rterra2023Matrix <- ez.matrify(rterra2023Filt, species.name = "ASVHeader",
                                  site.name = "sample", abundance = "asvAbsoluteAbundance")
    
    hist(colSums(rterra2023Matrix), breaks = 50)
    
    # Create a matrix where the replicates for sites are combined
    
    rterra2023LCSite <- rterra2023Filt %>%
        dplyr::select(LCSite, ASVHeader, asvAbsoluteAbundance) %>%
        group_by(LCSite, ASVHeader) %>%
        summarise(abundance = sum(asvAbsoluteAbundance))
    
    rterraSite2022Matrix <- ez.matrify(rterra2023LCSite, species.name = "ASVHeader",
                                      site.name = "LCSite", abundance = "abundance")
    
    #test to make sure everything got in
    any((colSums(rterraSite2022Matrix)-colSums(rterra2023Matrix)) != 0 )
    
    
    
    