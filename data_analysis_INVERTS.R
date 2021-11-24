## 2. Data analysis::Biodiversity 

## This file takes the processed eDNA data and conducts the biodiversity data
## analysis. This file specifically addresses the vertebrates using the Riaz primers.

## There were two protocls followed, one using 2g of soil and the other using 20g.

## ------ Setup --------------------------------------------


source("data_processing_IMPLEMENT.R")
source("../../../RCode/R_Scripts/PlotTaxaKD.R") # this code is found in my github here: https://github.com/kdyson/R_Scripts
source("../../../RCode/R_Scripts/repeat_multipatt.R") # ditto
library(vegan)
library(indicspecies)
library(TITAN2)
library(lme4)
library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)
detach("package:raster")

rm(list = ls()[grepl(ls(), pattern = "riaz|gill|zea")])

## ----- Summary of eDNA data ------------------------------

## eDNA reads found; total number of species variants, etc. 
invertData$total_reads <- rowSums(invertData[, grepl(x = colnames(invertData), "SFX[0-9]")])
invertData_2$total_reads <- rowSums(invertData_2[, grepl(x = colnames(invertData_2), "SFX[0-9]")])
invertData_20$total_reads <- rowSums(invertData_20[, grepl(x = colnames(invertData_20), "SFX[0-9]")])

totalReadsInvert <- sum(invertData$total_reads)
totalReadsInvert_2 <- sum(invertData_2$total_reads)
totalReadsInvert_20 <- sum(invertData_20$total_reads)

totalMOTUInvert <- length(unique(invertData$id[invertData$total_reads>0]))
totalMOTUInvert_2 <- length(unique(invertData_2$id[invertData_2$total_reads>0]))
totalMOTUInvert_20 <- length(unique(invertData_20$id[invertData_20$total_reads>0]))

totalUniqueTaxon <- length(unique(invertData$family))
totalUniqueTaxon_2 <- length(unique(invertData_2$family))
totalUniqueTaxon_20 <- length(unique(invertData_20$family))

totalReadsInvert_2/totalReadsInvert

colSums(invertData_2[ , grepl(x = colnames(invertData_2), pattern = "-")])
sum(rowSums(invertData_2[ , grepl(x = colnames(invertData_2), pattern = "-")])>0)



## _______ Key species questions and indicators ________________

## ----- Key Species present -----------------------------------
# What keystone species are present in the samples?

# Cocoa key species
cocoaKeyInvert_2 <- invert_2_cocoaSiteSpecies[, colnames(invert_2_cocoaSiteSpecies) %in% keyInvertList]
cocoaKeyInvert_20 <- invert_20_cocoaSiteSpecies[, colnames(invert_20_cocoaSiteSpecies) %in% keyInvertList]

cocoaKeyInvert_2_Table <- tibble(
    id = colnames(cocoaKeyInvert_2),
    `Cocoa Mean Abundance` = cocoaKeyInvert_2 %>%
        lapply(., mean, MARGIN = 2) %>% 
        unlist() %>% round(4),
    `Cocoa Standard Deviation` = cocoaKeyInvert_2 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)

cocoaKeyInvert_20_Table <- tibble(
    id = colnames(cocoaKeyInvert_20),
    `Cocoa Mean Abundance` = cocoaKeyInvert_20 %>%
        lapply(., mean, MARGIN = 2) %>% 
        unlist() %>% round(4),
    `Cocoa Standard Deviation` = cocoaKeyInvert_20 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)


# Pasture key species
pastureKeyInvert_2 <- invert_2_pastureSiteSpecies[, colnames(invert_2_pastureSiteSpecies) %in% keyInvertList]
pastureKeyInvert_20 <- invert_20_pastureSiteSpecies[, colnames(invert_20_pastureSiteSpecies) %in% keyInvertList]

pastureKeyInvert_2_Table <-  tibble(
    id = colnames(pastureKeyInvert_2),
    `Pasture Mean Abundance` = pastureKeyInvert_2 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4),
    `Pasture Standard Deviation` = pastureKeyInvert_2 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)

pastureKeyInvert_20_Table <-  tibble(
    id = colnames(pastureKeyInvert_20),
    `Pasture Mean Abundance` = pastureKeyInvert_20 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4),
    `Pasture Standard Deviation` = pastureKeyInvert_20 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)

# Forests
forestKeyInvert_2 <-
    invert_2_siteSpecies[grepl(x = rownames(invert_2_siteSpecies), pattern = "SFX004-01F"), colnames(invert_2_siteSpecies) %in% keyInvertList]

forestKeyInvert_20 <-
    invert_20_siteSpecies[grepl(x = rownames(invert_20_siteSpecies), pattern = "SFX004-01F"), colnames(invert_20_siteSpecies) %in% keyInvertList]
forestKeyInvert_20 <- forestKeyInvert_20[ , colSums(forestKeyInvert_20) > 0]

forestKeyInvert_2_Table <-  tibble(
    id = colnames(forestKeyInvert_2),
    `Forest Mean Abundance` = forestKeyInvert_2 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4),
    `Forest Standard Deviation` = forestKeyInvert_2 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)
forestKeyInvert_2_Table <- forestKeyInvert_2_Table[ forestKeyInvert_2_Table$`Forest Mean Abundance` > 0 , ] 



forestKeyInvert_20_Table <-  tibble(
    id = colnames(forestKeyInvert_20),
    `Forest Mean Abundance` = forestKeyInvert_20 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4),
    `Forest Standard Deviation` = forestKeyInvert_20 %>%
        lapply(., sd) %>%
        unlist() %>% round(2)
)

# create a table to present the data


keyInvert_2_Table <- tibble(id = keyInvertList) %>%
    left_join(keyInvertSpecies[, c("id", "species")], "id") %>%
    left_join(cocoaKeyInvert_2_Table, "id") %>%
    left_join(pastureKeyInvert_2_Table, "id") %>%
    left_join(forestKeyInvert_2_Table, "id")



keyInvert_20_Table <- tibble(
    id = keyInvertList
) %>%
    left_join(keyInvertSpecies[ , c("id", "species")], "id") %>%
    left_join(cocoaKeyInvert_20_Table, "id") %>%
    left_join(pastureKeyInvert_20_Table, "id") %>%
    left_join(forestKeyInvert_20_Table, "id")


remove(cocoaKeyInvert_2_Table, cocoaKeyInvert_20_Table,
       pastureKeyInvert_2_Table, pastureKeyInvert_20_Table,
       forestKeyInvert_2_Table, forestKeyInvert_20_Table)

## ----- Given what we know, use 2g going forward --------------

# Key species tables:
cocoaKeySpecies <- cocoaKeyInvert_2
pastureKeySpecies <- pastureKeyInvert_2
forestKeySpecies <- forestKeyInvert_2

remove(cocoaKeyInvert_2, cocoaKeyInvert_20,
       pastureKeyInvert_2, pastureKeyInvert_20,
       forestKeyInvert_2, forestKeyInvert_20)




# Species x Site tables:
cocoaSpeciesSite <- invert_2_cocoa
pastureSpeciesSite <- invert_2_pasture
allSpeciesSite <- invert_2_ReadOnly


# invert_2_ReadOnly
# invert_20_ReadOnly
# invert_2_siteSpecies


# Site x Species table
cocoaSiteSpecies <- invert_2_cocoaSiteSpecies
pastureSiteSpecies <- invert_2_pastureSiteSpecies
allSiteSpecies <- invert_2_siteSpecies

allSiteSpecies$landUse <- ifelse(grepl("[0-9][0-9]C",rownames(allSiteSpecies)),
                                 yes = "Cocoa",
                                 no = ifelse(grepl("[0-9][0-9]P", rownames(allSiteSpecies)),
                                             yes = "Pasture",
                                             no = "Forest"))



rm(list = ls()[grepl(ls(), pattern = "invert_2|_20")])



## ----- Biodiversity Indicator 4 -------------------------------
## Number of keystone/priority species due to intervention.

# setup
cocoaKSPresent <-
    cocoaKeySpecies %>%
    select(which(colSums(.) > 0)) %>% colnames()
pastureKSPresent <-
    pastureKeySpecies %>%
    select(which(colSums(.) > 0)) %>% colnames()
forestKSPresent <-
    forestKeySpecies %>%
    select(which(colSums(.) > 0)) %>% colnames()
allKSPresent <- allSiteSpecies[, colnames(allSiteSpecies) %in% keyInvertList] %>%
    colnames()

## Which key species are in cocoa and not in pasture? 
keyInCocoaOnly <- cocoaKSPresent[!(cocoaKSPresent %in% pastureKSPresent)]

keyInCocoaOnly_Table <- left_join(tibble(id = keyInCocoaOnly),
                                  keyInvertSpecies)

## Which key species are in pasture and not in cocoa?
keyInPastureOnly <- pastureKSPresent[!(pastureKSPresent %in% cocoaKSPresent)]

keyInPastureOnly_Table <- left_join(tibble(id = keyInPastureOnly),
                                  keyInvertSpecies)

## Which species are in both?
keyInBoth <- cocoaKSPresent[cocoaKSPresent %in% pastureKSPresent]

## Which species are in neither?
keyInNeither <- keyInvertList[!(keyInvertList %in% c(cocoaKSPresent, pastureKSPresent))]

## Which species are in forests Only
keyInForestOnly <- forestKSPresent[forestKSPresent %in% keyInNeither]

keyInNone2g <- keyInNeither[!(keyInNeither %in% forestKSPresent)]

## ----- Indicator species analysis --------------------------- 
## Indicator species analysis: which species are particularly associated with
## intervention sites versus counterfactual sites?

# Doing this without the forest ones.



temp <- allSiteSpecies %>%
    #    filter(landUse != "Forest") %>%
    select(which(colnames(allSiteSpecies) %in% keyInvertList))

indicatorSpecies <- repeat.multipatt(
    matrix.name = temp,
    cluster.name = allSiteSpecies$landUse,#[allSiteSpecies$landUse != "Forest"],
    quiet = F,
    p.cutoff = 0.05
)[[1]]

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

cocoaKeySpecies$updatedPlot <- rownames(cocoaKeySpecies)
cocoaKeySpecies <-
    cocoaKeySpecies %>% left_join(select(siteLookup, c(updatedPlot, site, system)))

pastureKeySpecies$updatedPlot <- rownames(pastureKeySpecies)
pastureKeySpecies <-
    pastureKeySpecies %>% left_join(select(siteLookup, c(updatedPlot, site, system)))



tbSpeciesChange <- tibble(
    siteNames = c(cocoaKeySpecies$updatedPlot, pastureKeySpecies$updatedPlot),
    siteField = c(cocoaKeySpecies$site, pastureKeySpecies$site),
    siteType = c(cocoaKeySpecies$system, pastureKeySpecies$system),
    keySpeciesCount = c(
        specnumber(select(cocoaKeySpecies, where(is.numeric))),
        specnumber(select(pastureKeySpecies, where(is.numeric)))),
    keySpeciesAbundance = c(rowSums(select(cocoaKeySpecies, where(is.numeric))), rowSums(select(pastureKeySpecies, where(is.numeric))))
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
test_keySpeciesCountLMER <- car::Anova(keySpeciesCountLMER)


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
anova(keySpeciesAbundanceLMER)
lmerTest::rand(keySpeciesAbundanceLMER)
summary(keySpeciesAbundanceLMER)
test_keySpeciesAbundanceLMER <- car::Anova(keySpeciesAbundanceLMER)



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
cocoaSpeciesRichness <- specnumber(cocoaSpeciesSite[-1] , MARGIN = 2)

mean(cocoaSpeciesRichness)
sd(cocoaSpeciesRichness)

# Calculate Species Richness for pasture
pastureSpeciesRichness <- specnumber(pastureSpeciesSite[-1] , MARGIN = 2)

mean(pastureSpeciesRichness)
sd(pastureSpeciesRichness)

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
test_allSpeciesRichnessLMER <- car::Anova(allSpeciesRichnessLMER)
lmerTest::rand(allSpeciesRichnessLMER)
summary(allSpeciesRichnessLMER)


## Now let's make some graphs: first, species accumulation curves.

cocoaSpecAccum <- specaccum(cocoaSiteSpecies)
pastureSpecAccum <- specaccum(pastureSiteSpecies)

plot(
    pastureSpecAccum,
    ci.type = "poly",
    col = "red",
    lwd = 2,
    ci.lty = 0,
    ci.col = "pink",
    main = "Species Accumulation",
    ylab = "Expected (mean) species richness"
)

plot(
    cocoaSpecAccum,
    ci.type = "poly",
    col = "blue",
    lwd = 2,
    ci.lty = 0,
    ci.col = "lightblue",
    add = TRUE
)

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

allAbundanceLMER <- lmerTest::lmer(allAbundance ~
                                   siteType +
                                   (1 | siteField),
                               data = tbSpeciesChange,
                               REML = TRUE)
anova(allAbundanceLMER)

lmerTest::rand(allAbundanceLMER)
summary(allAbundanceLMER)
test_allSpeciesAbundanceLMER <- car::Anova(allAbundanceLMER)


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
    relativeDominance = absoluteDominance / sum(select(allSiteSpecies,!(landUse)))#,
    #lowAbundance = apply(select(allSiteSpecies,!(landUse)), MARGIN = 1, FUN = log_modulo_skewness)
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
test_shannonLMER <- car::Anova(shannonLMER)

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




# ESR

effectiveSRLMER <- lme4::lmer(effectiveSR ~
                                  siteType +
                                  (1 | siteField),
                              data = tbSpeciesChange,
                              REML = TRUE)
anova(effectiveSRLMER)
lmerTest::rand(effectiveSRLMER)
summary(effectiveSRLMER)
test_effectiveSRLMER<-car::Anova(effectiveSRLMER)

tbSpeciesChange %>%
    ggplot() +
    geom_boxplot(aes(siteField, effectiveSR),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, effectiveSR, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Effective Species Richness", x = "Field Site Code")




# ... and the beta diversity.

#Whittaker index
allWhittaker <- betadiver(select(allSiteSpecies, !(landUse)), method = "w")

# Bray-Curtis
allBC <- vegdist(select(allSiteSpecies, !(landUse)))



alphaDiversity <- mean(specnumber(select(allSiteSpecies, !(landUse))))
gammaDiversity <- length(specnumber(select(allSiteSpecies, !(landUse)), MARGIN = 2) > rare)
betaDiversity <- gammaDiversity/alphaDiversity

betadisper(d = allBC, group = allSiteSpecies$landUse)
plot(betadisper(d = allBC, group = allSiteSpecies$landUse), main = "Invertebrate PCoA by Field Type (Bray Curtis)", sub = "")

nmds <- vegan::metaMDS(dplyr::select(allSiteSpecies, !(landUse)), distance = "bray", trymax = 999, autotransform = T)
plot(nmds)

## ----- END OF CODE ---------------------------