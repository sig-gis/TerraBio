# This analyis is based on soil samples collected in June at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and determines three
# things: 
# 1. Is Longmire's buffer or desiccant a better preservative to use in
# the field, where "better" is more community data and coverage returned?
# 2. How many replicates are needed per sample? 
# 3. How much soil should be collected from each site?

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


## ----- Data cleaning and setup -------------------------------
# Remove columns we don't need for analysis.
head(hortaAllASV)
hortaAllASV <- hortaAllASV[ , lookupColnames$Keep == "Y"]
head(hortaAllASV)

# Remove ASVs that don't meet criteria.
hortaAllASV <- hortaAllASV[ hortaAllASV$primerExpectedLength == "in range", ]


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

remove(hortaAllASV, hortaString)

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

## ----- Vol Testing | Replicates ------------------------------

#Create a curve within each letter sample for the three replicates. Basically
#trying to figure out how many new ATVs does each replicate contribute?


## ----- Vol Testing | Total Volume ----------------------------


## ----- CODE END ----------------------------------------------