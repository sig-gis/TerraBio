# This analysis is based on soil samples collected in September 2022 at Cafe Apui
# and processed by EcoMol.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis (Cafe Apui TB indicators)

## ----- Data ingestion ---------------------------------------

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
source("../../../RCode/R_Scripts/triplet_fixer.R") # from my github repository
source("../../../RCode/R_Scripts/repeat_multipatt.R") # ditto

lookupColnames <- read.csv("CafeApui/lookupColnames.csv")
lookupSitenames <- read.csv("CafeApui/lookupSitenames.csv")

apuiAllASV <- read.csv("CafeApui/Apui_2022_bioinfo_results.csv",
                        stringsAsFactors = F,
                        col.names = lookupColnames$TB_ColName)
    apuiAllASV$ASVHeader <- str_sub(apuiAllASV$ASVHeader, 2, -1) # strips the > from the ASV header


# script variable definitions
minlibrarySize = 5000
minRelativeAbund = 0.05
minAbsoluteAbund = 5
readOrigin = "merged"

phylum = c("Arthropoda")
#phylum = c("Annelida", "Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca") #worms and insects
#phylum = c("Rhodophyta", "Streptophyta", "Chlorophyta", "Bacillariophyta") #plants and algae
#phylum = c("Rotifera", "Tubulinea", "Oomycota", "Gastrotricha", "Discosea", "Endomyxa", "Ascomycota") #microorganisms and fungi

#Reitmeier, S., Hitch, T.C., Treichel, N., Fikas, N., Hausmann, B., Ramer-Tait,
#A.E., Neuhaus, K., Berry, D., Haller, D., Lagkouvardos, I. and Clavel, T.,
#2021. Handling of spurious sequences affects the outcome of high-throughput 16S
#rRNA gene amplicon profiling. ISME Communications, 1(1), pp.1-12.



## ----- Data cleaning  -------------------------------
# Remove columns we don't need for analysis.
head(apuiAllASV)
apuiAllASV <- apuiAllASV[ , lookupColnames$Keep == "Y"]
head(apuiAllASV)

# Remove replicates with library size as in Rule of Thumb
# approaches described in: Cao, Q., Sun, X., Rajesh, K., Chalasani, N., Gelow,
# K., Katz, B., Shah, V.H., Sanyal, A.J. and Smirnova, E., 2021. Effects of rare
# microbiome taxa filtering on statistical analysis. Frontiers in microbiology,
# 11, p.607325.

removedSites <- unique(apuiAllASV$sample[apuiAllASV$sampleTotalAbd <= minlibrarySize])
apuiAllASV <- apuiAllASV[ !(apuiAllASV$sample %in% removedSites) , ]
    stopifnot(length(unique(apuiAllASV$sample[apuiAllASV$sampleTotalAbd <= minlibrarySize])) == 0)

    print(removedSites)
    remove(removedSites)

# Remove ASVs that don't meet criteria.
apuiAllASV <- apuiAllASV[ apuiAllASV$primerExpectedLength == "in range" & 
                              apuiAllASV$phylumBLASTn %in% phylum , ]

# Sum of each read origin reads
apuiSummary <- apuiAllASV %>% group_by(readOrigin) %>% 
    summarise(
        totalReads = sum(asvAbsoluteAbundance),
        numberRepl = length(unique(fileName)),
        avgReadsRepl = totalReads/numberRepl,
        uniqueASVs = length(unique(ASVHeader)),
        avgASVs = uniqueASVs/numberRepl,
        uniqueGenus = length(unique(genusBLASTn)),
        uniqueIdent = length(unique(identification)),
        avgPseudo = mean(blastPseudoScore)
    )

# Select Read Origin (merged, forward, or backward)
apuiAllASV <- apuiAllASV[ apuiAllASV$readOrigin == readOrigin , ]

# Remove ASVs below minimum relative abundance on sample
plot(
    apuiAllASV$asvAbsoluteAbundance,
    apuiAllASV$relativeAbundanceOnSample,
    ylim = c(0, .5),
    xlim = c(0, 50)
)
abline(h = minRelativeAbund, col = "red")
abline(v = minAbsoluteAbund, col = "red")

print("Number of ASV under min relative abundance:")
print(sum(apuiAllASV$relativeAbundanceOnSample < minRelativeAbund))
print("Number of ASV under min absolute abundance:")
print(sum(apuiAllASV$asvAbsoluteAbundance < minAbsoluteAbund))
print(sum(apuiAllASV$asvAbsoluteAbundance > minAbsoluteAbund))


apuiAllASV <- apuiAllASV[ apuiAllASV$asvAbsoluteAbundance > minAbsoluteAbund, ]

## ----- Create Matrices for analysis ------------------

# Changing replicate names to be more readable
apuiAllASV <- merge(apuiAllASV, lookupSitenames, by.x = "sample", by.y = "EM_Sample")

# Now we need to split out the comparison samples and the silica only samples (main analysis).
apuiPairedASV <- apuiAllASV [ apuiAllASV$Paired == "Y" ,  ]
apuiASV <- apuiAllASV[ apuiAllASV$Storage == "Silica" , ]
apuiBufferASV <- apuiAllASV[ apuiAllASV$Storage == "Buffer" ,]

# Now create Site Species matrices for each.
apuiPairedMatrix <- ez.matrify(apuiPairedASV, species.name = "ASVHeader", site.name = "sampleID", abundance = "asvAbsoluteAbundance")

apuiMatrix <- ez.matrify(apuiASV, species.name = "ASVHeader", site.name = "sampleID", abundance = "asvAbsoluteAbundance")

apuiBufferMatrix <- ez.matrify(apuiBufferASV, species.name = "ASVHeader", site.name = "sampleID", abundance = "asvAbsoluteAbundance")


# remove columns with 0s that might result from split
apuiMatrix <- apuiMatrix [ , colSums(apuiMatrix) > 0 ]
apuiPairedMatrix <- apuiPairedMatrix[ , colSums(apuiPairedMatrix) > 0 ]


## ----- Buffer vs. Silica for Apui data ------------------------
## A limited subset of the samples 

# plot all ASV absolute abundance by buffer/silica
ggplot(apuiPairedASV, aes(x = Storage, y = asvAbsoluteAbundance)) + 
    geom_violin() + geom_jitter()

# plot amount of DNA recovered from each sample by buffer/silica
ggplot(apuiPairedASV, aes(x = Storage, y = dnaConcentration)) + 
    geom_boxplot()

# plot 'quality' measure for each sample by buffer/silica
ggplot(apuiPairedASV, aes(x = Storage, y = dnaPurity)) + 
    geom_boxplot()





# plot total ASV absolute abundance grouped by sample for buffer/silica
volumeSampleASV <- apuiPairedASV %>% group_by(farmNumber, treatment, Storage) %>% 
    summarise(totalReads = sum(asvAbsoluteAbundance), countASV = length(unique(ASVHeader)))

ggplot(volumeSampleASV, aes(x = Storage, y = totalReads)) + 
    geom_boxplot() + geom_jitter()

# plot total number of ASVs grouped by sample for buffer/silica
ggplot(volumeSampleASV, aes(x = Storage, y = countASV)) + 
    geom_boxplot() + geom_jitter()


remove(volumeSampleASV)
