## 1. Data processing::Biodiversity
## This file takes eDNA data obtained from soil samples and prepares it for further analysis.

## Load packages

library(dplyr)


## groupings

allSites <- "^id$|best_match|G[0-9][0-9]|P[0-9][0-9]"
cocoaSites <- "^id$|best_match|G[0-9][0-9]"
pastureSites <- "^id$|best_match|P[0-9][0-9]"
keySpecies <- "LC049653|KC193777|NC_009592|AB120717|nonsense"
    # this should be an identifier for whatever keySpecies we end up using. This
    # is just a random assortment for now
keySpeciesList <- strsplit(keySpecies, "\\|")[[1]]
speciesIDColumns <- "id$|best_match|species_list"




## Set variables
rare <- 2 # how many times do you need to detect a species to count it?

## Import eDNA data

motuData <- read.csv("BNA4.All_MOTUs_.csv", sep = ";", quote = "\'")

readsOnly <- motuData[ , grepl(allSites, colnames(motuData))]

cocoaOnly <- motuData[ , grepl(cocoaSites, colnames(motuData))]

pastureOnly <- motuData[ , grepl(pastureSites, colnames(motuData))]

## Create species x site and site x species tables for analysis. different
## packages want them in different ways.

cocoaSpeciesSite <- 
    cocoaOnly[,-1] %>% 
    group_by(best_match) %>% 
    summarise(across(G01:G45, sum)) 

cocoaSpeciesSite <- cocoaSpeciesSite[ rowSums(cocoaSpeciesSite[-1]) > 0, ]

cocoaSiteSpecies <- as.data.frame(cocoaSpeciesSite, row.names = cocoaSpeciesSite$best_match)
cocoaSiteSpecies <- as.data.frame(t(cocoaSiteSpecies[-1]))
colnames(cocoaSiteSpecies) <- cocoaSpeciesSite$best_match


pastureSpeciesSite <- 
    pastureOnly[,-1] %>% 
    group_by(best_match) %>% 
    summarise(across(P01:P19, sum)) 

pastureSpeciesSite <- pastureSpeciesSite[ rowSums(pastureSpeciesSite[-1]) > 0, ]

pastureSiteSpecies <- as.data.frame(pastureSpeciesSite, row.names = pastureSpeciesSite$best_match)
pastureSiteSpecies <- as.data.frame(t(pastureSiteSpecies[-1]))
colnames(pastureSiteSpecies) <- pastureSpeciesSite$best_match

allSpeciesSite <- 
    readsOnly %>% 
    group_by(best_match) %>% 
    summarise(across(G01:P19, sum)) 

allSpeciesSite <- allSpeciesSite[ rowSums(allSpeciesSite[-1]) > 0, ]

allSiteSpecies <- as.data.frame(allSpeciesSite, row.names = allSpeciesSite$best_match)
allSiteSpecies <- as.data.frame(t(allSiteSpecies[-1]))
colnames(allSiteSpecies) <- allSpeciesSite$best_match
allSiteSpecies$landUse <- ifelse(grepl("G[0-9][0-9]",rownames(allSiteSpecies)), 
                                 yes = "Cocoa",
                                 no = ifelse(grepl("P[0-9][0-9]", rownames(allSiteSpecies)),
                                 yes = "Pasture",
                                 no = "Unknown"))
    


# Create a species lookup table
speciesLookup <- motuData[, grepl(speciesIDColumns, colnames(motuData))]
    


## Import site level data

# Buffered land use characteristics

# Data collected in the field

# External datasets