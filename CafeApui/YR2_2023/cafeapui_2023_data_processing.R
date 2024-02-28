# This analysis is based on soil samples collected at Cafe Apui and processed by
# EcoMol. There are four land uses: I = intervention (agroforestry), R =
# control/forest, CF = pasture, CF2 = juquira (successional regrowth). The 2023
# analysis includes YR1 2022 and YR 2 2023 data.

# This file take the eDNA species table provided by EcoMol and prepares them for
# analysis.

# Code written Feb. 2024 by Karen Dyson


## ----- Data ingestion ---------------------------------------

# libraries
library(ggplot2)
library(gridExtra)
library(tidyr)
library(vegan)
library(stringr)
library(dplyr)

# source functions
source("../../../RCode/R_Scripts/triplet_fixer.R")
source("functions.R")

# script variable definitions
minlibrarySize = 5000
minAbsoluteAbund = 5
rare = 50
readOrigin = "merged"

#Reitmeier, S., Hitch, T.C., Treichel, N., Fikas, N., Hausmann, B., Ramer-Tait,
#A.E., Neuhaus, K., Berry, D., Haller, D., Lagkouvardos, I. and Clavel, T.,
#2021. Handling of spurious sequences affects the outcome of high-throughput 16S
#rRNA gene amplicon profiling. ISME Communications, 1(1), pp.1-12.

phylum = c("Arthropoda")
#phylum = c("Annelida", "Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca") #worms and insects



## ingest data

## common data
lookupColnames <- read.csv("lookupColnames.csv")
lookupSitenames2022 <- read.csv("CafeApui/YR1_2022/lookupSitenames2022.csv", strip.white = T)
lookupSitenames2023 <- read.csv("CafeApui/YR2_2023/lookupSitenames2023.csv", strip.white = T)


# 2022 eDNA data import

temp <- lookupColnames[order(lookupColnames$EM_OrderCApui2022), ] 
temp <- temp$TB_ColName[which(!is.na(temp$EM_ColNameCApui2022) & !is.na(temp$EM_OrderCApui2022))]

apui2022Raw <- read.csv("CafeApui/YR1_2022/Apui_2022_bioinfo_results.csv",
                       stringsAsFactors = F,
                       col.names = temp)
apui2022Raw$ASVHeader <- str_sub(apui2022Raw$ASVHeader, 2, -1) # strips the > from the ASV header

# 2023 eDNA data import
temp <- lookupColnames[order(lookupColnames$EM_OrderCApui2023), ] 
temp <- temp$TB_ColName[which(!is.na(temp$EM_ColNameCApui2023) & !is.na(temp$EM_OrderCApui2023))]

apui2023Raw <- read.csv("CafeApui/YR2_2023/Apui_2023_bioinfo_results.csv",
                         stringsAsFactors = F,
                         col.names = temp)
apui2023Raw$ASVHeader <- str_sub(apui2023Raw$ASVHeader, 2, -1) # strips the > from the ASV header

remove(temp)

# quick diagnostics on 2022 and 2023 data import

ecoMolSummary(apui2022Raw)
ecoMolSummary(apui2023Raw)

apui2022Raw %>%
    group_by(phylumBLASTn) %>%
    summarise(sum(asvAbsoluteAbundance))

apui2023Raw %>%
    group_by(phylumBLASTn) %>%
    summarise(sum(asvAbsoluteAbundance))

apui2023Raw %>%
    group_by(identificationMaxTaxon) %>%
    summarise(sum(asvAbsoluteAbundance))


## ----- Data cleaning and setup | 2022 -------------------------------

# Remove columns we don't need for analysis.
head(apui2022Raw)
apui2022Filt <- apui2022Raw[ , lookupColnames$TB_ColName[which(lookupColnames$TB_ColName %in% colnames(apui2022Raw) &
                                                                                    lookupColnames$Keep == "Y")]]
head(apui2022Filt)

# Remove replicates with small library size as in Rule of Thumb
# approaches described in: Cao, Q., Sun, X., Rajesh, K., Chalasani, N., Gelow,
# K., Katz, B., Shah, V.H., Sanyal, A.J. and Smirnova, E., 2021. Effects of rare
# microbiome taxa filtering on statistical analysis. Frontiers in microbiology,
# 11, p.607325.

removedSites2022 <- unique(apui2022Filt$sample[apui2022Filt$sampleTotalAbd <= minlibrarySize])
apui2022Filt <- apui2022Filt[ !(apui2022Filt$sample %in% removedSites2022) , ]
stopifnot(length(unique(apui2022Filt$sample[apui2022Filt$sampleTotalAbd <= minlibrarySize])) == 0)

print(removedSites2022)

# Remove ASVs that don't meet length and phylum criteria.
apui2022Filt <- apui2022Filt[ apui2022Filt$primerExpectedLength == "in range" & 
                                  apui2022Filt$phylumBLASTn %in% phylum , ]

# Sum of each read origin reads
apuiSummary2022 <- apui2022Filt %>% group_by(readOrigin) %>% 
    summarise(
        totalReads = sum(asvAbsoluteAbundance),
        numberRepl = length(unique(sample)),
        avgReadsRepl = totalReads/numberRepl,
        uniqueASVs = length(unique(ASVHeader)),
        avgASVs = uniqueASVs/numberRepl,
        uniqueGenus = length(unique(genusBLASTn)),
        uniqueIdent = length(unique(identificationFiltered)),
        avgPseudo = mean(blastPseudoScore)
    )

# Select Read Origin (merged, forward, or backward)
apui2022Filt <- apui2022Filt[ apui2022Filt$readOrigin == readOrigin , ]

# add information about the treatments to the ASV tables

apui2022Filt <- lookupSitenames2022 %>%
    select(EM_Sample, farmName, treatment, storage, replSampleID, shortSampleID) %>%
    right_join(., y = apui2022Filt, 
               by = c("EM_Sample" = "sample")
               )


## ----- Data cleaning and setup | 2023 -------------------------------

# Fix typo in the sample column, where EM135c5_CF3 should be EM135c5_8CF3

apui2023Raw$sample[ apui2023Raw$sample == "EM135c5_CF3"] <- "EM135c5_8CF3"


# Remove columns we don't need for analysis.
head(apui2023Raw)
apui2023Filt <- apui2023Raw[ , lookupColnames$TB_ColName[which(lookupColnames$TB_ColName %in% colnames(apui2023Raw) &
                                                                   lookupColnames$Keep == "Y")]]
head(apui2023Filt)

# Remove replicates with small library size as in Rule of Thumb
# approaches described in: Cao, Q., Sun, X., Rajesh, K., Chalasani, N., Gelow,
# K., Katz, B., Shah, V.H., Sanyal, A.J. and Smirnova, E., 2021. Effects of rare
# microbiome taxa filtering on statistical analysis. Frontiers in microbiology,
# 11, p.607325.

removedSites2023 <- unique(apui2023Filt$sample[apui2023Filt$sampleTotalAbd <= minlibrarySize])
apui2023Filt <- apui2023Filt[ !(apui2023Filt$sample %in% removedSites2023) , ]
stopifnot(length(unique(apui2023Filt$sample[apui2023Filt$sampleTotalAbd <= minlibrarySize])) == 0)

print(removedSites2023)

# Remove ASVs that don't meet length and phylum criteria.
apui2023Filt <- apui2023Filt[ apui2023Filt$primerExpectedLength == "in range" & 
                                  apui2023Filt$phylumBLASTn %in% phylum , ]

# Sum of each read origin reads (only merged provided in 2023 data)
apuiSummary2023 <- apui2023Filt %>% group_by(readOrigin) %>% 
    summarise(
        totalReads = sum(asvAbsoluteAbundance),
        numberRepl = length(unique(sample)),
        avgReadsRepl = totalReads/numberRepl,
        uniqueASVs = length(unique(ASVHeader)),
        avgASVs = uniqueASVs/numberRepl,
        uniqueGenus = length(unique(genusBLASTn)),
        uniqueIdent = length(unique(identificationFiltered)),
        avgPseudo = mean(blastPseudoScore)
    )


# add information about the treatments to the ASV tables

apui2023Filt <- lookupSitenames2023 %>%
    select(EM_Sample, farmName, treatment, storage, replSampleID, shortSampleID) %>%
    right_join(., y = apui2023Filt, 
               by = c("EM_Sample" = "sample")
    )



## ------ Identify common sequences | 2022 & 2023 ------------------------------------------

# This identifies sequences in both years' data--this suggests it is probably a
# "real" sequence and not a fluke. It also facilitates the joint PCAs
commonASV <- tibble(
    ASVSequence = apui2023Raw$ASVSequence[apui2023Raw$ASVSequence %in% apui2022Raw$ASVSequence]
)

commonASV <- merge(x = commonASV,
                   y = apui2022Raw[ , colnames(apui2022Raw) %in% c("ASVSequence", "ASVHeader")],
                   by = "ASVSequence",
                   all.x = T
) %>% unique()
colnames(commonASV)[2] <- "ASVHeader2022"
commonASV <- merge(x = commonASV,
                   y = apui2023Raw[ , colnames(apui2023Raw) %in% c("ASVSequence", "ASVHeader")],
                   by = "ASVSequence",
                   all.x = T
) %>% unique()
colnames(commonASV)[3] <- "ASVHeader2023"
commonASV$lookup <- paste0(commonASV$ASVHeader2022, "=",commonASV$ASVHeader2023)


# ----- Data quality checks | 2022 & 2023 -------------------------------

# Compare 2022 and 2023 DNA concentration.

p0 <- ggplot(subset(lookupSitenames2022, storage == "Silica"),
             aes(x = treatment, y = dnaConcentration)) + 
    geom_boxplot() + geom_point(color = "green") +
    ylim(5,80) + xlab("LC type (silica 2022)")

p1 <- ggplot(subset(lookupSitenames2022, storage == "Buffer"),
             aes(x = treatment, y = dnaConcentration)) + 
    geom_boxplot() + geom_point(color = "darkgreen") +
    ylim(5,80) + xlab("LC type (buffer 2022)")

p2 <- ggplot(lookupSitenames2023,
             aes(x = treatment, y = dnaConcentration)) + 
    geom_boxplot() + geom_point() +
    ylim(5,80) + xlab("LC type (silica 2023)")


grid.arrange(p0, p1, p2, nrow=2)

# Compare 2022 and 2023 DNA purity.

p1 <- ggplot(subset(lookupSitenames2022, storage == "Silica"),
             aes(x = treatment, y = dnaPurity)) + 
    geom_boxplot() + geom_point(color = "darkgreen") +
    ylim(1.3,2.2) + xlab("LC type (silica 2022)")

p2 <- ggplot(lookupSitenames2023,
             aes(x = treatment, y = dnaPurity)) + 
    geom_boxplot() + geom_point() +
    ylim(1.3,2.2) + xlab("LC type (silica 2023)")


grid.arrange(p1, p2, ncol=2)    

remove(p0, p1, p2)

# Look at ASV absolute abundance

ggplot(apui2022Filt, aes(x = treatment, y = asvAbsoluteAbundance)) + 
    geom_boxplot(outlier.colour = "red") +
    scale_y_log10()
#+ geom_jitter()

ggplot(apui2023Filt, aes(x = treatment, y = asvAbsoluteAbundance)) + 
    geom_boxplot(outlier.colour = "red") +
    scale_y_log10()
#+ geom_jitter()

# A few very high abundance ASVs 

#what are the weird ASVs?
apui2022Filt[apui2022Filt$asvAbsoluteAbundance > 7000, ]
apui2023Filt[apui2023Filt$asvAbsoluteAbundance > 7000, ]


# look at sample total abundance
ggplot(apui2022Filt, aes(x = treatment, y = sampleTotalAbd)) + 
    geom_boxplot() + geom_point()

ggplot(apui2023Filt, aes(x = treatment, y = sampleTotalAbd)) + 
    geom_boxplot() + geom_point() 


# ----- Filter for rare species | 2022 & 2023 -------------------------

# We want to remove rare species that could be errors 

# Remove ASVs below minimum relative abundance on sample. 
plot(
    apui2022Filt$asvAbsoluteAbundance,
    apui2022Filt$relativeAbundanceOnSample,
    ylim = c(0, .05),
    xlim = c(0, 50)
)
abline(v = minAbsoluteAbund, col = "red")

print("2022 filt. ASV < min absolute abundance:")
print(table(apui2022Filt$asvAbsoluteAbundance < minAbsoluteAbund))

apui2022Filt <- apui2022Filt[ apui2022Filt$asvAbsoluteAbundance >= minAbsoluteAbund, ]

# 2023: Remove ASVs below minimum relative abundance on sample. 
plot(
    apui2023Filt$asvAbsoluteAbundance,
    apui2023Filt$relativeAbundanceOnSample,
    ylim = c(0, .003),
    xlim = c(0, 50)
)
abline(v = minAbsoluteAbund, col = "red")

# 2023 has very low relative abundance on sample compared with apui 2022 and Cafe Apui 2022

print("2023 filt. ASV < min absolute abundance:")
print(table(apui2023Filt$asvAbsoluteAbundance < minAbsoluteAbund))

apui2023Filt <- apui2023Filt[ apui2023Filt$asvAbsoluteAbundance >= minAbsoluteAbund, ]




# ----- Create Matrices | 2022 ---------------------------------------

# Create a matrix with replicates as individual "sites"

apuiMatrix2022 <- ez.matrify(apui2022Filt, species.name = "ASVHeader",
                              site.name = "replSampleID", abundance = "asvAbsoluteAbundance")

hist(colSums(apuiMatrix2022), breaks = 50)
# Create a matrix where the replicates for land use types are combined

apuiType2022 <- apui2022Filt %>%
    dplyr::select(shortSampleID, ASVHeader, asvAbsoluteAbundance) %>%
    group_by(shortSampleID, ASVHeader) %>%
    summarise(abundance = sum(asvAbsoluteAbundance))

apuiMatrixType2022 <- ez.matrify(apuiType2022, species.name = "ASVHeader",
                                  site.name = "shortSampleID", abundance = "abundance")

#test to make sure everything got in
any((colSums(apuiMatrixType2022)-colSums(apuiMatrix2022)) != 0 )

# ----- Create Matrices | 2023 ---------------------------------------

# Create a matrix with replicates as individual "sites"

apuiMatrix2023 <- ez.matrify(apui2023Filt, species.name = "ASVHeader",
                              site.name = "replSampleID", abundance = "asvAbsoluteAbundance")

hist(colSums(apuiMatrix2023), breaks = 100)





# Create a matrix where the replicates for land use types are combined

apuiType2023 <- apui2023Filt %>%
    dplyr::select(shortSampleID, ASVHeader, asvAbsoluteAbundance) %>%
    group_by(shortSampleID, ASVHeader) %>%
    summarise(abundance = sum(asvAbsoluteAbundance))

apuiMatrixType2023 <- ez.matrify(apuiType2023, species.name = "ASVHeader",
                                  site.name = "shortSampleID", abundance = "abundance")

#test to make sure everything got in
any((colSums(apuiMatrixType2023)-colSums(apuiMatrix2023)) != 0 )


