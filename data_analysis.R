## 2. Data analysis::Biodiversity 

## This file takes the processed eDNA data and conducts the biodiversity data
## analysis (may need to split up for length)

## ------ Setup --------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#     BiocManager::install("microbiome")
library(microbiome)


source("data_processing.R")
source("../../../RCode/R_Scripts/PlotTaxaKD.R") # this code is found in my github here: https://github.com/kdyson/R_Scripts
source("../../../RCode/R_Scripts/repeat_multipatt.R") # ditto
library(dplyr)
library(vegan)
library(indicspecies)
library(TITAN2)
library(lme4)
library(ggplot2)


## ----- Summary of eDNA data ------------------------------

## eDNA reads found; total number of species variants, etc. 
totalReads <- sum(readsOnly[-c(1:2)])
totalMOTU <- length(unique(readsOnly$id))
totalUniqueBest <- length(unique(motuData$best_match))
totalUniqueTaxon <- length(unique(motuData$taxid))

## _______ Key species questions and indicators ________________

## ----- Key Species present -----------------------------------
# What keystone species are present in the samples?

# Cocoa key species
cocoaKeySpecies <- cocoaSiteSpecies[ , grepl(keySpecies,
                                             colnames(cocoaSiteSpecies))]

cocoaKSTable <- tibble(
    best_match = colnames(cocoaKeySpecies),
    `Cocoa Mean Abundance` = cocoaKeySpecies %>%
        lapply(., mean, MARGIN = 2) %>% 
        unlist() %>% round(4),
    `Cocoa Standard Deviation` = cocoaKeySpecies %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)


# Pasture key species
pastureKeySpecies <- pastureSiteSpecies[ , grepl(keySpecies,
                                                 colnames(pastureSiteSpecies))]

pastureKSTable <-  tibble(
    best_match = colnames(pastureKeySpecies),
    `Pasture Mean Abundance` = pastureKeySpecies %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4),
    `Pasture Standard Deviation` = pastureKeySpecies %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)

# create a table to present the data

keySpeciesTable <- tibble(
    best_match = keySpeciesList
) %>%
    left_join(cocoaKSTable, "best_match") %>%
    left_join(pastureKSTable, "best_match")# %>%
#    left_join(speciesLookup, "best_match")

# there are lots of superfulous matches, may want to change this later.

## ----- Biodiversity Indicator 4 -------------------------------
## Number of keystone/priority species due to intervention.

# setup
cocoaKSPresent <-
    select(cocoaKeySpecies,!starts_with("siteT")) %>%
    select(which(colSums(.) > 0)) %>% colnames()
pastureKSPresent <-
    select(pastureKeySpecies,!starts_with("siteT")) %>% 
    select(which(colSums(.) > 0)) %>% colnames()
allKSPresent <- allSiteSpecies[ , grepl(keySpecies,
                                              colnames(allSiteSpecies))] %>%
    colnames()

## Which key species are in cocoa and not in pasture? 
keyInCocoaOnly <- cocoaKSPresent[!(cocoaKSPresent %in% pastureKSPresent)]

## Which key species are in pasture and not in cocoa?
keyInPastureOnly <- pastureKSPresent[!(pastureKSPresent %in% cocoaKSPresent)]


## Which species are in both?
keyInBoth <- cocoaKSPresent[cocoaKSPresent %in% pastureKSPresent]
    
## Which species are in neither?
keyInNeither <- keySpeciesList[!(keySpeciesList %in% c(cocoaKSPresent, pastureKSPresent))]

## ----- Indicator species analysis --------------------------- 
## Indicator species analysis: which species are particularly associated with
## intervention sites versus counterfactual sites?

temp <- select(allSiteSpecies, !(landUse))

indicatorSpecies <- repeat.multipatt(matrix.name = temp,
                                     cluster.name = allSiteSpecies$landUse, 
                                     quiet = F,
                                     p.cutoff = 0.05)[[1]]

remove(temp)

indicatorSpecies <-
    indicatorSpecies[order(indicatorSpecies$groupname, indicatorSpecies$mean.stat, decreasing = TRUE), ] %>%
    select(groupname, species, frequency.sp, mean.stat, min.p.val, max.p.val, )







## ----- Question 1.8A ------------------------------------------
## Does the number of keystone/priority species (i.e. key pollinator species)
## change over time compared to counterfactual?

# Since we have a repeated measures design we need to use GLMM.

specnumber(cocoaKeySpecies)

# First create a number of keystone/priority species for each 

tbSpeciesChange <- tibble(
    siteNames = c(rownames(cocoaKeySpecies), rownames(pastureKeySpecies)),
    siteField = c(cocoaField, pastureField),
    siteType = c(rep("Cocoa", length(
        rownames(cocoaKeySpecies)
    )), rep("Pasture", length(
        rownames(pastureKeySpecies)
    ))),
    keySpeciesCount = c(
        specnumber(cocoaKeySpecies) + rnorm(length(rownames(cocoaKeySpecies))),
        specnumber(pastureKeySpecies) + rnorm(length(rownames(pastureKeySpecies)))
    ),
    keySpeciesAbundance = c(rowSums(cocoaKeySpecies), rowSums(pastureKeySpecies))
)

tbSpeciesChange$siteField <- factor(tbSpeciesChange$siteField, levels = unique(tbSpeciesChange$siteField))

temp1 <- tbSpeciesChange %>% group_by(siteField) %>% summarise(meanCount = mean(keySpeciesCount))
temp2 <- tbSpeciesChange %>% group_by(siteField) %>% summarise(meanAbundance = mean(keySpeciesAbundance))
temp2$x <- seq(0.75,nrow(temp2) + 0.25, by = 1)
temp2$xend <- seq(1.25,nrow(temp2) + 0.25,1)



tbSpeciesChange <- left_join(tbSpeciesChange, temp1) %>% left_join(temp2)
tbSpeciesChange$siteField <- factor(tbSpeciesChange$siteField, levels = temp1$siteField)

remove(temp1, temp2)


library(lmerTest) # this may mask rand.


keySpeciesCountLMER <- lme4::lmer(keySpeciesCount ~
                                      siteType +
                                      (1 | siteField),
                                  data = tbSpeciesChange,
                                  REML = TRUE)
anova(keySpeciesCountLMER)
lmerTest::rand(keySpeciesCountLMER)
summary(keySpeciesCountLMER)
car::Anova(keySpeciesCountLMER)

sjPlot::plot_model(keySpeciesCountLMER, show.values=T, show.p = T)


# Create graph where x-axis has each field, with individual points for each
# plot. y-axis should be number of species. Add in means for fields, means & 95%
# CI for cocoa/pasture, and possibly significance as well.


plotTBSpeciesCount <- tbSpeciesChange %>%
    ggplot() +
    geom_boxplot(aes(siteField, keySpeciesCount),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, keySpeciesCount, color = siteType),
        width = 0.1,
        height = 0
    ) +
    # geom_segment(aes(
    #     x = x,
    #     y = meanCount,
    #     xend = xend,
    #     yend = meanCount
    # ),
    # color = "black",
    # size = 1.25) +
    #
    # facet_grid(cols=vars(siteType), space = "free", scales = "free_x") +
    labs(color = "Field Type", y = "Key Species Richness", x = "Field Site Code")




## ----- Question 1.8B -------------------------------------------
## Does the abundance of keystone/priority species (i.e. key pollinator species)
## change over time compared to counterfactual? & NP Test for Biodiversity
## Indicator 5 (ABF-KPI-8): Change in abundance of keystone/ priority species
## due to interventions.

# Since we have a repeated measures design we need to use GLMM.


keySpeciesAbundanceLMER <- lme4::lmer(keySpeciesAbundance ~
                                      siteType +
                                      (1 | siteField),
                                  data = tbSpeciesChange,
                                  REML = TRUE)
anova(keySpeciesCountLMER)
lmerTest::rand(keySpeciesCountLMER)
summary(keySpeciesCountLMER)
car::Anova(keySpeciesCountLMER)




# Create graph where x-axis has each field, with individual points for each
# plot. y-axis should be abundance of species. Add in means for fields, means & 95%
# CI for cocoa/pasture, and possibly significance as well.


plotTBSpeciesAbundance <- tbSpeciesChange %>%
    ggplot() +
    geom_boxplot(aes(siteField, keySpeciesAbundance),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, keySpeciesAbundance, color = siteType),
        width = 0.1,
        height = 0
    ) +
    # geom_segment(aes(
    #     x = x,
    #     y = meanAbundance,
    #     xend = xend,
    #     yend = meanAbundance
    # ),
    # color = "black",
    # size = 1.25) +
    
    #    facet_grid(cols=vars(siteType), space = "free", scales = "free_x") +
    labs(color = "Field Type", y = "Key Species Abundance", x = "Field Site Code")


#https://cmdlinetips.com/2019/05/how-to-highlight-select-data-points-with-ggplot2-in-r/
#https://ggplot2.tidyverse.org/reference/geom_segment.html

## _____ General species richness, abundance, & diversity questions ____________________

## ----- Question 1.5 --------------------------------------------

## Does species richness change over time
## compared to counterfactual? & Biodiversity Indicator 6: Change in species
## richness due to interventions. H0: Increase in species richness in cocoa
## fields compared with pasture. H1: No change or decrease in species richness.

# note this differs from the above because it is all species richness, not just key species richness.

# Calculate species richness for cocoa and pasture, then add it to the existing tbSpeciesChange table

    # Calculate Species Richness for cocoa
    cocoaSpeciesRichness <- colSums(cocoaSpeciesSite[-1] >= rare)
    mean(cocoaSpeciesRichness)
    sd(cocoaSpeciesRichness)
    cocoaTotalSR <- sum(rowSums(cocoaSpeciesSite[-1]) >= rare)
    # total number of species is 183 (total detected 192). Average per site is 68.3 SD 9.20
    
    # Calculate Species Richness for pasture
    pastureSpeciesRichness <- colSums(pastureSpeciesSite[-1] >= rare)
    mean(pastureSpeciesRichness)
    sd(pastureSpeciesRichness)
    pastureTotalSR <- sum(rowSums(pastureSpeciesSite[-1]) >= rare) # if rare = 0 then it is equal to nrow
    # total number of species is 155 out of 163. Average per site is 65.7 SD 9.71
    
    temp <- tibble(siteNames = c(labels(cocoaSpeciesRichness), labels(pastureSpeciesRichness)), 
                   allSpeciesRichness = c(cocoaSpeciesRichness, pastureSpeciesRichness)
    )
    
    tbSpeciesChange <- left_join(tbSpeciesChange, temp)
    remove(temp)

# Since we have a repeated measures design we need to use GLMM.


allSpeciesRichnessLMER <- lme4::lmer(allSpeciesRichness ~
                                          siteType +
                                          (1 | siteField),
                                      data = tbSpeciesChange,
                                      REML = TRUE)
anova(allSpeciesRichnessLMER)
lmerTest::rand(allSpeciesRichnessLMER)
summary(allSpeciesRichnessLMER)
car::Anova(allSpeciesRichnessLMER)


## Now let's make some graphs: first, species accumulation curves.

cocoaSpecAccum <- specaccum(cocoaSiteSpecies)
pastureSpecAccum <- specaccum(pastureSiteSpecies)

plot(cocoaSpecAccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", 
     main = "Cocoa Species Accumulation")

plot(pastureSpecAccum, ci.type="poly", col="red", lwd=2, ci.lty=0, ci.col="pink", add=TRUE)

# now plot by field


tbSpeciesAllCount <- tbSpeciesChange %>%
    ggplot() +
    geom_boxplot(aes(siteField, allSpeciesRichness),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, allSpeciesRichness, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "All Species Richness", x = "Field Site Code")



## ----- Question 1.6 ------------------------------
# NP Test for Question 1.6: Does relative abundance of species change over time
# compared to counterfactual? H0: Increase in relative abundance of species on
# cocoa fields compared with pastures. H1: No change or decrease in relative
# abundance of species.

# Calculate species richness for cocoa and pasture, then add it to the existing tbSpeciesChange table

    # Calcuate relative abundance for cocoa
    cocoaAbundance <- rowSums(cocoaSiteSpecies)
    
            mean(cocoaAbundance)
            sd(cocoaAbundance)

    # Calculate relative abundance for pasture
    pastureAbundance <- rowSums(pastureSiteSpecies)
            
        mean(pastureAbundance)
        sd(pastureAbundance)

    temp <- tibble(siteNames = c(labels(cocoaAbundance), labels(pastureAbundance)), 
                   allAbundance = c(cocoaAbundance, pastureAbundance)
    )
    
    tbSpeciesChange <- left_join(tbSpeciesChange, temp)
    remove(temp)
        
        
# Test to compare

    allAbundanceLMER <- lme4::lmer(allAbundance ~
                                             siteType +
                                             (1 | siteField),
                                         data = tbSpeciesChange,
                                         REML = TRUE)
    anova(allAbundanceLMER)
    lmerTest::rand(allAbundanceLMER)
    summary(allAbundanceLMER)
    car::Anova(allAbundanceLMER)
    
    
# for specific species instead of totals; please add columns with that species' data.
    
    
    tbSpeciesAllAbundance <- tbSpeciesChange %>%
        ggplot() +
        geom_boxplot(aes(siteField, allAbundance),
                     outlier.shape = NA) +
        geom_jitter(aes(siteField, allAbundance, color = siteType),
                    width = 0.1,
                    height = 0) +
        labs(color = "Field Type", y = "All Species Abundance (reads)", x = "Field Site Code")
    




## ----- Question 1.7 --------------------------------------------
# Are there changes in Shannon’s diversity index (and others) over time compared
# to counterfactual? & Biodiversity Indicator 7: Change in biodiversity indices
# due to interventions.

# Here, we'll want to calculate the site-specific diversity metrics...

diversityMetrics <- tibble(
    siteNames = rownames(allSiteSpecies),
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
    
    tbSpeciesChange <- left_join(tbSpeciesChange, diversityMetrics)


# Test to compare site diversities between land use types
    # Test to compare
    # shannons.
    
    shannonLMER <- lme4::lmer(shannon ~
                                  siteType +
                                  (1 | siteField),
                              data = tbSpeciesChange,
                              REML = TRUE)
    anova(shannonLMER)
    lmerTest::rand(shannonLMER)
    summary(shannonLMER)
    car::Anova(shannonLMER)
    
    tbSpeciesChange %>%
        ggplot() +
        geom_boxplot(aes(siteField, shannon),
                     outlier.shape = NA) +
        geom_jitter(
            aes(siteField, shannon, color = siteType),
            width = 0.1,
            height = 0
        ) +
        labs(color = "Field Type", y = "Shannon Diversity", x = "Field Site Code")
    
    
    
    
    # pielou
    
    pielouLMER <- lme4::lmer(pielou ~
                                 siteType +
                                 (1 | siteField),
                             data = tbSpeciesChange,
                             REML = TRUE)
    anova(pielouLMER)
    lmerTest::rand(pielouLMER)
    summary(pielouLMER)
    car::Anova(pielouLMER)
    
    tbSpeciesChange %>%
        ggplot() +
        geom_boxplot(aes(siteField, pielou),
                     outlier.shape = NA) +
        geom_jitter(
            aes(siteField, pielou, color = siteType),
            width = 0.1,
            height = 0
        ) +
        labs(color = "Field Type", y = "Pielou", x = "Field Site Code")
    



# ... and the beta diversity.

#Whittaker index
allWhittaker <- betadiver(select(allSiteSpecies, !(landUse)), method = "w")

# Bray-Curtis
allBC <- vegdist(select(allSiteSpecies, !(landUse)))



alphaDiversity <- mean(specnumber(allSiteSpecies))
gammaDiversity <- length(specnumber(allSiteSpecies, MARGIN = 2) > rare)
betaDiversity <- gammaDiversity/alphaDiversity

betadisper(d = allBC, group = allSiteSpecies$landUse)
plot(betadisper(d = allBC, group = allSiteSpecies$landUse))


# note this doesnt account properly for the repeated measures.
adonis(allBC ~ allSiteSpecies$landUse)




## _____ Habitat framework ___________________________

# what species do we care about???

# Need environmental variables--NDVI or NDFI, secondary forest data, tree height, and elevation along with land use should be enough.
# Possible packages:
# https://onlinelibrary.wiley.com/doi/10.1111/ecog.02671
# 
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.01881
# 
# https://rspatial.org/sdm/SDM.pdf
# 
#SSDM: https://mran.microsoft.com/web/packages/SSDM/vignettes/SSDM.html
## ----- END OF CODE ---------------------------