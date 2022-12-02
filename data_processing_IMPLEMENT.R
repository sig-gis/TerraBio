## 1. Data processing::Biodiversity This file takes- eDNA data obtained from
## soil samples and prepares it for further analysis. Note that there were three
## primers used, one for mammals and two for invertebrates.

## Load packages

library(plyr)
library(dplyr)
library(vegan)

## Set variables
rare <- 1 # minimum times do you need to detect a species to include.
minAbun <- 50 # minimum number of reads for a species to include.
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

gilletData <- read.csv("PilotFarmData/Gillet_SFX.csv", sep = ",", quote = "\'")

riazData <- read.csv("PilotFarmData/Riaz_SFX.csv", sep = ",", quote = "\'")

zealeData <- read.csv("PilotFarmData/Zeale_SFX.csv", sep = ",", quote = "\'")


sum(rowSums(riazData[ , grepl(x = colnames(riazData), pattern = "_2g")])>0)
    # riaz has a whole bunch of extra 0s, for some reason.

sum(rowSums(gilletData[ , grepl(x = colnames(gilletData), pattern = "_2")])>0)
sum(rowSums(zealeData[ , grepl(x = colnames(zealeData), pattern = "_2")])>0)



## ----- Remove suspicious plots ----------------------------

riazData <- dplyr::select(riazData, !(SFX6_2g | SFX6_20g))
    # plot had an abnormally high amount of cow reads
sum(rowSums(riazData[ , grepl(x = colnames(riazData), pattern = "_2")])>0)
sum(rowSums(riazData[ , grepl(x = colnames(riazData), pattern = "_2g")])>0)

gilletData <- dplyr::select(gilletData, !(SFX39_2))
    # There is no SFX39

## ----- Remove bad phylums & 0s --------------------------

gilletData <- gilletData[gilletData$phylum == "Arthropoda", ]
zealeData <- zealeData[zealeData$phylum == "Arthropoda", ]

riazData <- riazData[riazData$total_reads > 0, ]
gilletData <- gilletData[gilletData$total_reads > 0 , ]
zealeData <- zealeData[zealeData$total_reads > 0 , ]


invertData <- invertBind(gilletData,zealeData)
sum(rowSums(invertData[ , grepl(x = colnames(invertData), pattern = "_")]) > 0)


# ----- Import key species tables -----------------------------------

keyVertSpecies <- read.csv("PilotFarmData/KeyVertebrates.csv")
keyInvertSpecies <- read.csv("PilotFarmData/KeyInvertebrates.csv")


# this should be an identifier for keySpecies
keyVertList <- keyVertSpecies$id[keyVertSpecies$endangered %in% 
                                     c("vulnerable", "decreasing", "extinct")]
vertList <- keyVertSpecies$id[keyVertSpecies$keySpecies == "yes"]
keyInvertList <- keyInvertSpecies$id[keyInvertSpecies$keySpecies == "yes"]


# Import lookup table, which has field codes and the codes the Salford lab used
siteLookup <- read.csv("PilotFarmData/Lookup_SFX.csv")

# ----- Split the OTU datasets in to 2 gram and 20 gram methods --------------

riazData_2 <-
    riazData[, grepl(x = colnames(riazData), pattern = "^id$|_2g|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]
riazData_20 <-
    riazData[, grepl(x = colnames(riazData), pattern = "^id$|_20g|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]

sum(rowSums(riazData[ , grepl(x = colnames(riazData), pattern = "_2g")])>0)
sum(rowSums(riazData_2[ , grepl(x = colnames(riazData_2), pattern = "_2g")])>0)
sum(rowSums(riazData_20[ , grepl(x = colnames(riazData_20), pattern = "_20g")])>0)

# update total abundances.
riazData_2$total_reads <- rowSums(riazData_2[ , grepl(x = colnames(riazData_2), pattern = "_2g")])
riazData_20$total_reads <- rowSums(riazData_20[ , grepl(x = colnames(riazData_20), pattern = "_20g")])




gilletData_2 <-
    gilletData[, grepl(x = colnames(gilletData), pattern = "^id$|_2$|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]
gilletData_20 <-
    gilletData[, grepl(x = colnames(gilletData), pattern = "^id$|_20|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]

sum(rowSums(gilletData[ , grepl(x = colnames(gilletData), pattern = "_2$")])>0)
sum(rowSums(gilletData_2[ , grepl(x = colnames(gilletData_2), pattern = "_2$")])>0)
sum(rowSums(gilletData_20[ , grepl(x = colnames(gilletData_20), pattern = "_20")])>0)

# update total abundances.
gilletData_2$total_reads <- rowSums(gilletData_2[ , grepl(x = colnames(gilletData_2), pattern = "_2")])
gilletData_20$total_reads <- rowSums(gilletData_20[ , grepl(x = colnames(gilletData_20), pattern = "_20")])





zealeData_2 <-
    zealeData[, grepl(x = colnames(zealeData), pattern = "^id$|_2$|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]
zealeData_20 <-
    zealeData[, grepl(x = colnames(zealeData), pattern = "^id$|_20|king|phy|cla|ord|fam|gen|spe|pid|eva|tot|seq")]

sum(rowSums(zealeData[ , grepl(x = colnames(zealeData), pattern = "_2$")])>0)
sum(rowSums(zealeData_2[ , grepl(x = colnames(zealeData_2), pattern = "_2$")])>0)
sum(rowSums(zealeData_20[ , grepl(x = colnames(zealeData_20), pattern = "_20")])>0)


# update total abundances.
zealeData_2$total_reads <- rowSums(zealeData_2[ , grepl(x = colnames(zealeData_2), pattern = "_2")])
zealeData_20$total_reads <- rowSums(zealeData_20[ , grepl(x = colnames(zealeData_20), pattern = "_20")])







## Replace the sample IDs used by the lab with the field ones.


colnames(riazData_2) <- mapvalues(colnames(riazData_2),
                                  from = siteLookup$riazCode2,
                                  to = siteLookup$updatedPlot)
colnames(riazData_20) <- mapvalues(colnames(riazData_20),
                                   from = siteLookup$riazCode20,
                                   to = siteLookup$updatedPlot)



colnames(gilletData_2) <- mapvalues(colnames(gilletData_2),
                                     from = siteLookup$labCode2,
                                     to = siteLookup$updatedPlot)

colnames(gilletData_20) <- mapvalues(colnames(gilletData_20),
                                    from = siteLookup$labCode20,
                                    to = siteLookup$updatedPlot)

colnames(zealeData_2) <- mapvalues(colnames(zealeData_2),
                                    from = siteLookup$labCode2,
                                    to = siteLookup$updatedPlot)
colnames(zealeData_20) <- mapvalues(colnames(zealeData_20),
                                    from = siteLookup$labCode20,
                                    to = siteLookup$updatedPlot)







## ----- Implement rare and minAbun filters ------------------------------

## Now that protocols have been separated, remove from each MOTUs with fewer
## reads than the chosen minimum abundance.


riazData_2 <- riazData_2[riazData_2$total_reads > minAbun, ]
riazData_20 <- riazData_20[riazData_20$total_reads > minAbun, ]


gilletData_2 <- gilletData_2[gilletData_2$total_reads > minAbun , ]
gilletData_20 <- gilletData_20[gilletData_20$total_reads > minAbun , ]

zealeData_2 <- zealeData_2[zealeData_2$total_reads > minAbun , ]
zealeData_20 <- zealeData_20[zealeData_20$total_reads > minAbun , ]



## And now remove MOTUs with fewer observations than minimum.

riazData_2 <- riazData_2[rowSums(riazData_2[, grepl(x = colnames(riazData_2), pattern = "-")] > 0) >= rare, ]
riazData_20 <- riazData_20[rowSums(riazData_20[, grepl(x = colnames(riazData_20), pattern = "-")] > 0) >= rare, ]


gilletData_2 <- gilletData_2[rowSums(gilletData_2[, grepl(x = colnames(gilletData_2), pattern = "-")] > 0) >= rare, ]
gilletData_20 <- gilletData_20[rowSums(gilletData_20[, grepl(x = colnames(gilletData_20), pattern = "-")] > 0) >= rare, ]

zealeData_2 <- zealeData_2[rowSums(zealeData_2[, grepl(x = colnames(zealeData_2), pattern = "-")] > 0) >= rare, ]
zealeData_20 <- zealeData_20[rowSums(zealeData_20[, grepl(x = colnames(zealeData_20), pattern = "-")] > 0) >= rare, ]


## Now smush the two invert datasets together (Gillet and Zeale primers)
any(gilletData$sequence %in% zealeData$sequence)

invertData_2 <- invertBind(gilletData_2,zealeData_2)
invertData_20 <- invertBind(gilletData_20,zealeData_20)

sum(rowSums(invertData_2[ , grepl(x = colnames(invertData_2), pattern = "-")])>0)

sum(rowSums(gilletData_20[ , grepl(x = colnames(gilletData_20), pattern = "-")])>0)
sum(rowSums(zealeData_20[ , grepl(x = colnames(zealeData_20), pattern = "-")])>0)
sum(rowSums(invertData_20[ , grepl(x = colnames(invertData_20), pattern = "-")])>0)

# good, matches as expected

if(cleanup == T) {
    remove(gilletData,
           gilletData_2,
           gilletData_20,
           zealeData,
           zealeData_2,
           zealeData_20)
}


## ------- Subset the eDNA data -------------------------------------



# Reads only tables for site/species matrices
riazData_2_ReadsOnly <- riazData_2[ , grepl(allSites, colnames(riazData_2))]
riazData_20_ReadsOnly <- riazData_20[ , grepl(allSites, colnames(riazData_20))]


invert_2_ReadsOnly <- invertData_2[ , grepl(allSites, colnames(invertData_2))]
invert_20_ReadsOnly <- invertData_20[ , grepl(allSites, colnames(invertData_20))]


# Cocoa fields
riazData_2_cocoaOnly <- riazData_2[ , grepl(cocoaSites, colnames(riazData_2))]
riazData_20_cocoaOnly <- riazData_20[ , grepl(cocoaSites, colnames(riazData_20))]


invert_2_cocoaOnly <- invertData_2[ , grepl(cocoaSites, colnames(invertData_2))]
invert_20_cocoaOnly <- invertData_20[ , grepl(cocoaSites, colnames(invertData_20))]

# Pastures
riazData_2_pastureOnly <- riazData_2[ , grepl(pastureSites, colnames(riazData_2))]
riazData_20_pastureOnly <- riazData_20[ , grepl(pastureSites, colnames(riazData_20))]

invert_2_pastureOnly <- invertData_2[ , grepl(pastureSites, colnames(invertData_2))]
invert_20_pastureOnly <- invertData_20[ , grepl(pastureSites, colnames(invertData_20))]



# forests
riazData_2_forestOnly <- riazData_2[ , grepl(forestSites, colnames(riazData_2))]
riazData_20_forestOnly <- riazData_20[ , grepl(forestSites, colnames(riazData_20))]

invert_2_forestOnly <- invertData_2[ , grepl(forestSites, colnames(invertData_2))]
invert_20_forestOnly <- invertData_20[ , grepl(forestSites, colnames(invertData_20))]






## ----- Species x Site and Site x Species ---------------------------------

## Create species x site and site x species tables for analysis. different
## packages want them in different ways.

# speciesSite <- function(inputOTU, rare = 0) {
#     temp <- inputOTU[rowSums(inputOTU[-1]) > rare,]
#     return(temp)
# }

siteSpecies <- function(inputOTU) {
    
    temp <- as.data.frame(inputOTU, row.names = inputOTU$id)
    temp <- as.data.frame(t(temp[-1]))
    colnames(temp) <- inputOTU$id
    
    return(temp)
    
}


# Inverts

    invert_2_siteSpecies <- siteSpecies(invert_2_ReadsOnly)
    invert_20_siteSpecies <- siteSpecies(invert_20_ReadsOnly)
    

    sum(rowSums(invert_2_ReadsOnly[ , grepl(x = colnames(invert_2_ReadsOnly), pattern = "-")])>0)
    sum(colSums(invert_2_siteSpecies[grepl(x = rownames(invert_2_siteSpecies), pattern = "-") , ])>0)
    
    
    
# Cocoa fields

    invert_2_cocoaSiteSpecies <- siteSpecies(invert_2_cocoaOnly)
    invert_20_cocoaSiteSpecies <- siteSpecies(invert_20_cocoaOnly)

    
# Pastures

    invert_2_pastureSiteSpecies <- siteSpecies(invert_2_pastureOnly)
    invert_20_pastureSiteSpecies <- siteSpecies(invert_20_pastureOnly)



# Vertebrates
    
    riazData_2_siteSpecies <- siteSpecies(riazData_2_ReadsOnly)
    riazData_20_siteSpecies <- siteSpecies(riazData_20_ReadsOnly)
    
    riazData_2_cocoaSiteSpecies <- siteSpecies(riazData_2_cocoaOnly)
    riazData_20_cocoaSiteSpecies <- siteSpecies(riazData_20_cocoaOnly)
    
    riazData_2_pastureSiteSpecies <- siteSpecies(riazData_2_pastureOnly)
    riazData_20_pastureSiteSpecies <- siteSpecies(riazData_20_pastureOnly)
    

    length(unique(riazData_2_ReadsOnly$id[rowSums(riazData_2_ReadsOnly[-1])>0]))



## ----- MISC ----------------------------
# field plot groups.

cocoaField <- (siteLookup$site[siteLookup$system == "COCOA"])

pastureField <- (siteLookup$site[siteLookup$system == "PASTURE"])

forestField <- (siteLookup$site[siteLookup$system == "FOREST"])



# allSiteSpecies$landUse <- ifelse(grepl("G[0-9][0-9]",rownames(allSiteSpecies)), 
#                                  yes = "Cocoa",
#                                  no = ifelse(grepl("P[0-9][0-9]", rownames(allSiteSpecies)),
#                                              yes = "Pasture",
#                                              no = "Forest"))
# 


## Import site level data.
# 
# library(raster)
# 
# # Need lat and long for species distribution models.
# sampledPlots <- shapefile("PilotFarmData/sampledPlots.shp", verbose = T)
# 
# plotNames <- unique(sampledPlots$ID)
# 
# sampledPlotsTable <- tibble(plotName = sampledPlots$ID,
#                             xCentroid = sampledPlots$xCentroid,
#                             yCentroid = sampledPlots$yCentroid)
# 
# detach("package:raster")


## few metrics for quality checks
invertData %>% select(SFX10_2:SFX9_20) %>% colSums()
invertData %>% select(SFX10_2:SFX9_20) %>% sum()
invertData_2 %>% select(`SFX237-04P-03`:`SFX237-04P-02`) %>% sum()

riazData %>% select(SFX10_2g:SFX9_20g) %>% colSums()
riazData %>% select(SFX10_2g:SFX9_20g) %>% sum()
riazData_2 %>% select(`SFX237-04P-03`:`SFX237-04P-02`) %>% sum()


## END
