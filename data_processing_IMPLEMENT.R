## 1. Data processing::Biodiversity
## This file takes eDNA data obtained from soil samples and prepares it for further analysis. Note that there were three primers used, one for mammals and two for invertebrates. 

## Load packages

library(plyr)
library(dplyr)
library(raster)
library(vegan)

## Set variables
rare <- 0 # how many times do you need to detect a species to count it?
perm <- 99999 # number of permutations to use
cleanup <- TRUE

## define groupings

allSites <- "^id$|SFX[0-9][0-9]"
cocoaSites <- "^id$|[0-9][0-9]C"
pastureSites <- "^id$|[0-9][0-9]P"
forestSites <-  "^id$|[0-9][0-9]F"


invertBind <- function(inputOTU1, inputOTU2) {
    
    # fill in non-overlapping columns with 0s
    inputOTU1[setdiff(names(inputOTU2), names(inputOTU1))] <- 0
    inputOTU2[setdiff(names(inputOTU1), names(inputOTU2))] <- 0
    temp <- rbind(inputOTU1, inputOTU2)
    return(temp)
    
}


# ----- Import the motu datasets -----------------------------------

gilletData <- read.csv("FarmData/Gillet_SFX.csv", sep = ",", quote = "\'")

riazData <- read.csv("FarmData/Riaz_SFX.csv", sep = ",", quote = "\'")

zealeData <- read.csv("FarmData/Zeale_SFX.csv", sep = ",", quote = "\'")


sum(rowSums(riazData[ , grepl(x = colnames(riazData), pattern = "_2g")])>0)


sum(rowSums(gilletData[ , grepl(x = colnames(gilletData), pattern = "_")])>0)
sum(rowSums(zealeData[ , grepl(x = colnames(zealeData), pattern = "_")])>0)



## ----- Remove suspicious plot ----------------------------

riazData <- dplyr::select(riazData, !(SFX6_2g))
sum(rowSums(riazData[ , grepl(x = colnames(riazData), pattern = "_2")])>0)
sum(rowSums(riazData[ , grepl(x = colnames(riazData), pattern = "_2g")])>0)

gilletData <- gilletData[ , colnames(gilletData) != "SFX39_2"]


## ----- Remove bad phylums & 0s ---------------------------------

gilletData <- gilletData[gilletData$phylum == "Arthropoda", ]
zealeData <- zealeData[zealeData$phylum == "Arthropoda", ]

gilletData <- gilletData[gilletData$total_reads > 0 , ]
zealeData <- zealeData[zealeData$total_reads > 0 , ]


invertData <- invertBind(gilletData,zealeData)
sum(rowSums(invertData[ , grepl(x = colnames(invertData), pattern = "_")])>0)


# ----- Import key species tables -----------------------------------
keyVertSpecies <- read.csv("FarmData/KeyVertebrates.csv")
keyInvertSpecies <- read.csv("FarmData/KeyInvertebrates.csv")


# this should be an identifier for keySpecies
keyVertList <- keyVertSpecies$id[keyVertSpecies$keySpecies == "yes"]
keyInvertList <- keyInvertSpecies$id[keyInvertSpecies$keySpecies == "yes"]


# Import lookup table, which has field codes and the codes the US lab used
siteLookup <- read.csv("FarmData/Lookup_SFX.csv")

# ----- Split the OTU datasets in to 2 gram and 20 gram methods --------------

riazData_2 <-
    riazData[, grepl(x = colnames(riazData), pattern = "^id$|_2g|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]
riazData_20 <-
    riazData[, grepl(x = colnames(riazData), pattern = "^id$|_20g|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]

sum(rowSums(riazData_2[ , grepl(x = colnames(riazData_2), pattern = "_2g")])>0)
sum(rowSums(riazData_20[ , grepl(x = colnames(riazData_20), pattern = "_20g")])>0)


gilletData_2 <-
    gilletData[, grepl(x = colnames(gilletData), pattern = "^id$|_2$|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]
gilletData_20 <-
    gilletData[, grepl(x = colnames(gilletData), pattern = "^id$|_20|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]

sum(rowSums(gilletData[ , grepl(x = colnames(gilletData), pattern = "_2$")])>0)
sum(rowSums(gilletData_2[ , grepl(x = colnames(gilletData_2), pattern = "_2$")])>0)
sum(rowSums(gilletData_20[ , grepl(x = colnames(gilletData_20), pattern = "_20")])>0)


zealeData_2 <-
    zealeData[, grepl(x = colnames(zealeData), pattern = "^id$|_2$|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]
zealeData_20 <-
    zealeData[, grepl(x = colnames(zealeData), pattern = "^id$|_20|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]

sum(rowSums(zealeData[ , grepl(x = colnames(zealeData), pattern = "_2$")])>0)
sum(rowSums(zealeData_2[ , grepl(x = colnames(zealeData_2), pattern = "_2$")])>0)
sum(rowSums(zealeData_20[ , grepl(x = colnames(zealeData_20), pattern = "_20")])>0)


## Replace the sample IDs used by the lab with the field ones.

colnames(gilletData_2) <- mapvalues(colnames(gilletData_2),
                                     from = siteLookup$labCode2,
                                     to = siteLookup$updatedPlot)

colnames(gilletData_20) <- mapvalues(colnames(gilletData_20),
                                    from = siteLookup$labCode20,
                                    to = siteLookup$updatedPlot)

colnames(riazData_2) <- mapvalues(colnames(riazData_2),
                                    from = siteLookup$riazCode2,
                                    to = siteLookup$updatedPlot)
colnames(riazData_20) <- mapvalues(colnames(riazData_20),
                                    from = siteLookup$riazCode20,
                                    to = siteLookup$updatedPlot)

colnames(zealeData_2) <- mapvalues(colnames(zealeData_2),
                                    from = siteLookup$labCode2,
                                    to = siteLookup$updatedPlot)
colnames(zealeData_20) <- mapvalues(colnames(zealeData_20),
                                    from = siteLookup$labCode20,
                                    to = siteLookup$updatedPlot)


invertData_2 <- invertBind(gilletData_2,zealeData_2)
invertData_20 <- invertBind(gilletData_20,zealeData_20)

sum(rowSums(invertData_2[ , grepl(x = colnames(invertData_2), pattern = "-")])>0)

sum(rowSums(gilletData_20[ , grepl(x = colnames(gilletData_20), pattern = "-")])>0)
sum(rowSums(zealeData_20[ , grepl(x = colnames(zealeData_20), pattern = "-")])>0)
sum(rowSums(invertData_20[ , grepl(x = colnames(invertData_20), pattern = "-")])>0)

# good, matches as expected


## ------- Subset the eDNA data -------------------------------------

any(gilletData$sequence %in% zealeData$sequence)

gilletData_2_ReadsOnly <- gilletData_2[ , grepl(allSites, colnames(gilletData_2))]
gilletData_20_ReadsOnly <- gilletData_20[ , grepl(allSites, colnames(gilletData_20))]

riazData_2_ReadsOnly <- riazData_2[ , grepl(allSites, colnames(riazData_2))]
riazData_20_ReadsOnly <- riazData_20[ , grepl(allSites, colnames(riazData_20))]

length(specnumber(riazData_2_ReadsOnly, MARGIN = 2) > rare)
length(specnumber(riazData_20_ReadsOnly, MARGIN = 2) > rare)


zealeData_2_ReadsOnly <- zealeData_2[ , grepl(allSites, colnames(zealeData_2))]
zealeData_20_ReadsOnly <- zealeData_20[ , grepl(allSites, colnames(zealeData_20))]

# Cocoa fields
gilletData_2_cocoaOnly <- gilletData_2[ , grepl(cocoaSites, colnames(gilletData_2))]
gilletData_20_cocoaOnly <- gilletData_20[ , grepl(cocoaSites, colnames(gilletData_20))]

riazData_2_cocoaOnly <- riazData_2[ , grepl(cocoaSites, colnames(riazData_2))]
riazData_20_cocoaOnly <- riazData_20[ , grepl(cocoaSites, colnames(riazData_20))]

zealeData_2_cocoaOnly <- zealeData_2[ , grepl(cocoaSites, colnames(zealeData_2))]
zealeData_20_cocoaOnly <- zealeData_20[ , grepl(cocoaSites, colnames(zealeData_20))]

# Pastures
gilletData_2_pastureOnly <- gilletData_2[ , grepl(pastureSites, colnames(gilletData_2))]
gilletData_20_pastureOnly <- gilletData_20[ , grepl(pastureSites, colnames(gilletData_20))]

riazData_2_pastureOnly <- riazData_2[ , grepl(pastureSites, colnames(riazData_2))]
riazData_20_pastureOnly <- riazData_20[ , grepl(pastureSites, colnames(riazData_20))]

zealeData_2_pastureOnly <- zealeData_2[ , grepl(pastureSites, colnames(zealeData_2))]
zealeData_20_pastureOnly <- zealeData_20[ , grepl(pastureSites, colnames(zealeData_20))]

# forests
gilletData_2_forestOnly <- gilletData_2[ , grepl(forestSites, colnames(gilletData_2))]
gilletData_20_forestOnly <- gilletData_20[ , grepl(forestSites, colnames(gilletData_20))]

riazData_2_forestOnly <- riazData_2[ , grepl(forestSites, colnames(riazData_2))]
riazData_20_forestOnly <- riazData_20[ , grepl(forestSites, colnames(riazData_20))]

zealeData_2_forestOnly <- zealeData_2[ , grepl(forestSites, colnames(zealeData_2))]
zealeData_20_forestOnly <- zealeData_20[ , grepl(forestSites, colnames(zealeData_20))]




## ----- Species x Site and Site x Species ---------------------------------

## Create species x site and site x species tables for analysis. different
## packages want them in different ways.

speciesSite <- function(inputOTU) {
    temp <- inputOTU[rowSums(inputOTU[-1]) > 0,]
    return(temp)
}

siteSpecies <- function(inputOTU) {
    
    temp <- as.data.frame(inputOTU, row.names = inputOTU$id)
    temp <- as.data.frame(t(temp[-1]))
    colnames(temp) <- inputOTU$id
    
    return(temp)
    
}


# Inverts

    gilletData_2_ReadsOnly <- speciesSite(gilletData_2_ReadsOnly)
    gilletData_20_ReadsOnly <- speciesSite(gilletData_20_ReadsOnly)

    zealeData_2_ReadsOnly <- speciesSite(zealeData_2_ReadsOnly)
    zealeData_20_ReadsOnly <- speciesSite(zealeData_20_ReadsOnly)

    sum(rowSums(gilletData_2_ReadsOnly[ , grepl(x = colnames(gilletData_2_ReadsOnly), pattern = "-")])>0)


    ## Combine the invertebrates
    
    invert_2_ReadOnly <- invertBind(gilletData_2_ReadsOnly, zealeData_2_ReadsOnly)
    invert_20_ReadOnly <- invertBind(gilletData_20_ReadsOnly, zealeData_20_ReadsOnly)
    invert_2_siteSpecies <- siteSpecies(invert_2_ReadOnly)
    invert_20_siteSpecies <- siteSpecies(invert_20_ReadOnly)
    
if(cleanup == TRUE){remove(gilletData_2_ReadsOnly, 
                           gilletData_20_ReadsOnly, 
                           zealeData_2_ReadsOnly, 
                           zealeData_20_ReadsOnly)}
    
    sum(rowSums(invert_2_ReadOnly[ , grepl(x = colnames(invert_2_ReadOnly), pattern = "-")])>0)
    sum(colSums(invert_2_siteSpecies[grepl(x = rownames(invert_2_siteSpecies), pattern = "-") , ])>0)
    
    
    
# Cocoa fields
    gilletData_2_cocoaOnly <- speciesSite(gilletData_2_cocoaOnly)
    gilletData_20_cocoaOnly <- speciesSite(gilletData_20_cocoaOnly)

    zealeData_2_cocoaOnly <- speciesSite(zealeData_2_cocoaOnly)
    zealeData_20_cocoaOnly <- speciesSite(zealeData_20_cocoaOnly)

    ## Combine the invertebrates
    
    invert_2_cocoa  <- invertBind(gilletData_2_cocoaOnly, zealeData_2_cocoaOnly)
    invert_20_cocoa <- invertBind(gilletData_20_cocoaOnly, zealeData_20_cocoaOnly)
    invert_2_cocoaSiteSpecies <- siteSpecies(invert_2_cocoa)
    invert_20_cocoaSiteSpecies <- siteSpecies(invert_20_cocoa)
    
    if(cleanup == TRUE){remove(gilletData_2_cocoaOnly, 
                               gilletData_20_cocoaOnly, 
                               zealeData_2_cocoaOnly, 
                               zealeData_20_cocoaOnly)}
    
# Pastures
    gilletData_2_pastureOnly <- speciesSite(gilletData_2_pastureOnly)
    gilletData_20_pastureOnly <- speciesSite(gilletData_20_pastureOnly)

    zealeData_2_pastureOnly <- speciesSite(zealeData_2_pastureOnly)
    zealeData_20_pastureOnly <- speciesSite(zealeData_20_pastureOnly)

    ## Combine the invertebrates
    
    invert_2_pasture <- invertBind(gilletData_2_pastureOnly, zealeData_2_pastureOnly)
    invert_20_pasture <- invertBind(gilletData_20_pastureOnly, zealeData_20_pastureOnly)
    invert_2_pastureSiteSpecies <- siteSpecies(invert_2_pasture)
    invert_20_pastureSiteSpecies <- siteSpecies(invert_20_pasture)
    
    if(cleanup == TRUE){remove(gilletData_2_pastureOnly, 
                               gilletData_20_pastureOnly, 
                               zealeData_2_pastureOnly, 
                               zealeData_20_pastureOnly)}
    
    

# Vertebrates
    
    riazData_2_ReadsOnly <- speciesSite(riazData_2_ReadsOnly)
    riazData_2_siteSpecies <- siteSpecies(riazData_2_ReadsOnly)
    
    riazData_20_ReadsOnly <- speciesSite(riazData_20_ReadsOnly)
    riazData_20_siteSpecies <- siteSpecies(riazData_20_ReadsOnly)
    
    riazData_2_cocoaOnly <- speciesSite(riazData_2_cocoaOnly)
    riazData_2_cocoaSiteSpecies <- siteSpecies(riazData_2_cocoaOnly)
    
    riazData_20_cocoaOnly <- speciesSite(riazData_20_cocoaOnly)
    riazData_20_cocoaSiteSpecies <- siteSpecies(riazData_20_cocoaOnly)
    
    riazData_2_pastureOnly <- speciesSite(riazData_2_pastureOnly)
    riazData_2_pastureSiteSpecies <- siteSpecies(riazData_2_pastureOnly)
    
    riazData_20_pastureOnly <- speciesSite(riazData_20_pastureOnly)
    riazData_20_pastureSiteSpecies <- siteSpecies(riazData_20_pastureOnly)
    

    length(unique(riazData_2_ReadsOnly$id[rowSums(riazData_2_ReadsOnly[-1])>0]))

# allSiteSpecies$landUse <- ifelse(grepl("G[0-9][0-9]",rownames(allSiteSpecies)), 
#                                  yes = "Cocoa",
#                                  no = ifelse(grepl("P[0-9][0-9]", rownames(allSiteSpecies)),
#                                              yes = "Pasture",
#                                              no = "Forest"))
# 


## ----- MISC ----------------------------
# field plot groups.

cocoaField <- (siteLookup$site[siteLookup$system == "COCOA"])

pastureField <- (siteLookup$site[siteLookup$system == "PASTURE"])

forestField <- (siteLookup$site[siteLookup$system == "FOREST"])




## Import site level data.

# Need lat and long for species distribution models.
sampledPlots <- shapefile("FarmData/sampledPlots.shp", verbose = T)

plotNames <- unique(sampledPlots$ID)

sampledPlotsTable <- tibble(plotName = sampledPlots$ID,
                            xCentroid = sampledPlots$xCentroid,
                            yCentroid = sampledPlots$yCentroid)



## END
