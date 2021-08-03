## 2. Data analysis::Biodiversity 

## This file takes the processed eDNA data and conducts the biodiversity data
## analysis (may need to split up for length)

## ----- Setup ------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("microbiome")

source("data_processing.R")
source("../../../RCode/R_Scripts/PlotTaxaKD.R") # this code is found in my github here: https://github.com/kdyson/R_Scripts
library(dplyr)
library(vegan)
library(indicspecies)
library(TITAN2)
library(microbiome)

## eDNA reads found; total number of species variants, etc. 
totalReads <- sum(readsOnly[-c(1:2)])
totalMOTU <- length(unique(readsOnly$id))
totalUniqueBest <- length(unique(motuData$best_match))
totalUniqueTaxon <- length(unique(motuData$taxid))


## ----- Question 1.5 --------------------------------------------

## Does species richness change over time
## compared to counterfactual? & Biodiversity Indicator 6: Change in species
## richness due to interventions. H0: Increase in species richness in cocoa
## fields compared with pasture. H1: No change or decrease in species richness.

# Calculate Species Richness for cocoa
cocoaSpeciesRichness <- colSums(cocoaSpeciesSite[-1] >= rare)
    mean(cocoaSpeciesRichness)
    sd(cocoaSpeciesRichness)
cocoaTotalSR <- sum(rowSums(cocoaSpeciesSite[-1]) >= rare)
# total number of species is 183 (total detected 192). Average per site is 68.3 SD 9.20

cocoaSpecAccum <- specaccum(cocoaSiteSpecies)
plot(cocoaSpecAccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", 
     main = "Cocoa Species Accumulation")


# Calculate Species Richness for pasture
pastureSpeciesRichness <- colSums(pastureSpeciesSite[-1] >= rare)
    mean(pastureSpeciesRichness)
    sd(pastureSpeciesRichness)
pastureTotalSR <- sum(rowSums(pastureSpeciesSite[-1]) >= rare) # if rare = 0 then it is equal to nrow
# total number of species is 155 out of 163. Average per site is 65.7 SD 9.71

pastureSpecAccum <- specaccum(pastureSiteSpecies)
plot(pastureSpecAccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", 
     main = "Pasture Species Accumulation")

# Test to compare
allSpeciesRichness <-
    tibble(
        landUse = c(rep("Cocoa", length(
            cocoaSpeciesRichness
        )), rep(
            "Pasture", length(pastureSpeciesRichness)
        )),
        speciesRichness = c(cocoaSpeciesRichness,
                            pastureSpeciesRichness)
    )

adonis(speciesRichness ~ landUse, data = allSpeciesRichness)





## ----- Question 1.6 ------------------------------
# NP Test for Question 1.6: Does relative abundance of species change over time
# compared to counterfactual? H0: Increase in relative abundance of species on
# cocoa fields compared with pastures. H1: No change or decrease in relative
# abundance of species.

# Calcuate relative abundance for cocoa
    cocoaKeySpecies <- cocoaSiteSpecies[ , grepl(keySpecies,
                                              colnames(cocoaSiteSpecies))]
    
    cocoaKeySpecies$siteTotal <- rowSums(cocoaKeySpecies)
        mean(cocoaKeySpecies$siteTotal)
        sd(cocoaKeySpecies$siteTotal)
#cocoaTotalSR <- sum(cocoaKeySpecies$siteTotal)


# Calculate relative abundance for pasture
    pastureKeySpecies <- pastureSiteSpecies[ , grepl(keySpecies,
                                                 colnames(pastureSiteSpecies))]
    
    pastureKeySpecies$siteTotal <- rowSums(pastureKeySpecies)
    mean(pastureKeySpecies$siteTotal)
    sd(pastureKeySpecies$siteTotal)
    #pastureTotalSR <- sum(pastureKeySpecies$siteTotal)

# Test to compare
allKeyAbundance <-
    tibble(
        landUse = c(rep("Cocoa", length(
            cocoaKeySpecies$siteTotal
        )), rep(
            "Pasture", length(pastureKeySpecies$siteTotal)
        )),
        speciesAbundance = c(cocoaKeySpecies$siteTotal,
                            pastureKeySpecies$siteTotal)
    )

adonis(speciesAbundance ~ landUse, data = allKeyAbundance)
# for specific species instead of totals; please add columns with that species' data.

## ----- Indicator species analysis --------------------------- 
## Indicator species analysis: which species are particularly associated with
## intervention sites versus counterfactual sites?

indicatorSpecies <- multipatt(select(allSiteSpecies, !(landUse)),
                        cluster = as.factor(allSiteSpecies$landUse))
summary(indicatorSpecies)

## -- Question 1.7 --------------------------------------------
# Are there changes in Shannon’s diversity index (and others) over time compared
# to counterfactual? & Biodiversity Indicator 7: Change in biodiversity indices
# due to interventions.

# Here, we'll want to calculate the site-specific diversity metrics...

diversityMetrics <- tibble(
    siteName = rownames(allSiteSpecies),
    landUse = allSiteSpecies$landUse,
    speciesRichness = rowSums(select(allSiteSpecies,!(landUse)) >= rare),
    shannon = vegan::diversity(
        select(allSiteSpecies,!(landUse)),
        index = "shannon",
        MARGIN = 1
    ),
    simpson =  vegan::diversity(
        select(allSiteSpecies,!(landUse)),
        index = "simpson",
        MARGIN = 1
    ),
    pielou = shannon / log(specnumber(select(
        allSiteSpecies,!(landUse)
    ))),
    effectiveSR = exp(shannon),
    absoluteDominance = apply(select(allSiteSpecies,!(landUse)), MARGIN = 1, max),
    relativeDominance = absoluteDominance / sum(select(allSiteSpecies,!(landUse))),
    lowAbundance = apply(select(allSiteSpecies,!(landUse)), MARGIN = 1, FUN = log_modulo_skewness)
)


# Test to compare site diversities between land use types
adonis(shannon ~ landUse, data = diversityMetrics)
    # shannons.




# ... and the beta diversity.

#Whittaker index
allWhittaker <- betadiver(select(allSiteSpecies, !(landUse)), method = "w")
# Bray-Curtis
allBC <- vegdist(select(allSiteSpecies, !(landUse)))

betadisper(d = allBC, group = allSiteSpecies$landUse)
plot(betadisper(d = allBC, group = allSiteSpecies$landUse))

adonis(allBC ~ allSiteSpecies$landUse)
