## 2. Data analysis::Biodiversity 

## This file takes the processed eDNA data and conducts the biodiversity data
## analysis. This file specifically addresses the invertebrates using the Gillet
## and Zeale primers.

## There were two protocols followed, one using 2g of soil and the other using 20g.

## Updated for pilot manuscript 5/27.

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
library(compositions)
library(zCompositions)
library(reshape2)
library(dplyr)


rm(list = ls()[grepl(ls(), pattern = "riaz|gill|zea")])

## ----- Summary of eDNA data ------------------------------

## note invertData etc. are without filtering... 2g and 20g protocols have been filtered.

## eDNA reads found; total number of species variants, etc. 
invertData$total_reads <- rowSums(invertData[, grepl(x = colnames(invertData), "SFX[0-9]")])
# Uncomment below if you don't update these in data processing.
# invertData_2$total_reads <- rowSums(invertData_2[, grepl(x = colnames(invertData_2), "SFX[0-9]")])
# invertData_20$total_reads <- rowSums(invertData_20[, grepl(x = colnames(invertData_20), "SFX[0-9]")])

totalReadsInvert_unfiltered <- sum(invertData$total_reads)
totalReadsInvert_filtered <- sum(invertData$total_reads[invertData$total_reads >= minAbun])
totalReadsInvert_2 <- sum(invertData_2$total_reads)
totalReadsInvert_20 <- sum(invertData_20$total_reads)

totalMOTUInvert_unfiltered <- length(unique(invertData$id[invertData$total_reads>0])) # need this check on total reads if e.g. minAbun = 0
totalMOTUInvert_filtered <- length(unique(invertData$id[invertData$total_reads >= minAbun]))

totalMOTUInvert_2 <- length(unique(invertData_2$id))
totalMOTUInvert_20 <- length(unique(invertData_20$id))

totalUniqueTaxon_unfiltered <- unique(invertData$family) %>% grep("Class_|Order_", ., invert = TRUE) %>% length()
totalUniqueTaxon_filtered <- unique(invertData$family[invertData$total_reads >= minAbun]) %>% grep("Class_|Order_", ., invert = TRUE) %>% length()

totalUniqueTaxon_2 <- unique(invertData_2$family) %>% grep("Class_|Order_", ., invert = TRUE) %>% length()
totalUniqueTaxon_20 <- unique(invertData_20$family) %>% grep("Class_|Order_", ., invert = TRUE) %>% length()






totalReadsInvert_2/totalReadsInvert_filtered

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
        unlist() %>% round(4) %>% unname(),
    `Cocoa Site Count` = cocoaKeyInvert_2 %>%
        lapply(., specnumber) %>% 
        unlist() %>% round(2) %>% unname()
)

cocoaKeyInvert_20_Table <- tibble(
    id = colnames(cocoaKeyInvert_20),
    `Cocoa Mean Abundance` = cocoaKeyInvert_20 %>%
        lapply(., mean, MARGIN = 2) %>% 
        unlist() %>% round(4) %>% unname(),
    `Cocoa Site Count` = cocoaKeyInvert_20 %>%
        lapply(., specnumber) %>% 
        unlist() %>% round(2) %>% unname()
)


# Pasture key species
pastureKeyInvert_2 <- invert_2_pastureSiteSpecies[, colnames(invert_2_pastureSiteSpecies) %in% keyInvertList]
pastureKeyInvert_20 <- invert_20_pastureSiteSpecies[, colnames(invert_20_pastureSiteSpecies) %in% keyInvertList]

pastureKeyInvert_2_Table <-  tibble(
    id = colnames(pastureKeyInvert_2),
    `Pasture Mean Abundance` = pastureKeyInvert_2 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4) %>% unname(),
    `Pasture Site Count` = pastureKeyInvert_2 %>%
        lapply(., specnumber) %>% 
        unlist() %>% round(2) %>% unname()
)

pastureKeyInvert_20_Table <-  tibble(
    id = colnames(pastureKeyInvert_20),
    `Pasture Mean Abundance` = pastureKeyInvert_20 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4) %>% unname(),
    `Pasture Site Count` = pastureKeyInvert_20 %>%
        lapply(., specnumber) %>% 
        unlist() %>% round(2) %>% unname()
)

# Forests
forestKeyInvert_2 <-
    invert_2_siteSpecies[grepl(x = rownames(invert_2_siteSpecies), pattern = "SFX004-01F"), colnames(invert_2_siteSpecies) %in% keyInvertList]




forestKeyInvert_2_Table <-  tibble(
    id = colnames(forestKeyInvert_2),
    `Forest Mean Abundance` = forestKeyInvert_2 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4) %>% unname(),
    `Forest Site Count` = forestKeyInvert_2 %>%
        lapply(., specnumber) %>% 
        unlist() %>% round(2) %>% unname()
)
forestKeyInvert_2_Table <- forestKeyInvert_2_Table[ forestKeyInvert_2_Table$`Forest Mean Abundance` > 0 , ] 



forestKeyInvert_20 <-
    invert_20_siteSpecies[grepl(x = rownames(invert_20_siteSpecies), pattern = "SFX004-01F"), colnames(invert_20_siteSpecies) %in% keyInvertList]
forestKeyInvert_20 <-
    forestKeyInvert_20[, colSums(forestKeyInvert_20) > 0]

forestKeyInvert_20_Table <-  tibble(
    id = colnames(forestKeyInvert_20),
    `Forest Mean Abundance` = forestKeyInvert_20 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4) %>% unname(),
    `Forest Site Count` = forestKeyInvert_20 %>%
        lapply(., specnumber) %>%
        unlist() %>% round(2) %>% unname()
)

# create a table to present the data

keyInvert_2_Table <- tibble(id = keyInvertList) %>%
    left_join(keyInvertSpecies[, c("id", "order", "family", "species")], "id") %>%
    left_join(cocoaKeyInvert_2_Table, "id") %>%
    left_join(pastureKeyInvert_2_Table, "id") %>%
    left_join(forestKeyInvert_2_Table, "id")

keyInvert_2_Table <- keyInvert_2_Table[ !(is.na(keyInvert_2_Table[ , 5])), ]

# keyInvert_20_Table <- tibble(
#     id = keyInvertList
# ) %>%
#     left_join(keyInvertSpecies[ , c("id", "order, "family", "species")], "id") %>%
#     left_join(cocoaKeyInvert_20_Table, "id") %>%
#     left_join(pastureKeyInvert_20_Table, "id") %>%
#     left_join(forestKeyInvert_20_Table, "id")


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
cocoaSpeciesSite <- invert_2_cocoaOnly
pastureSpeciesSite <- invert_2_pastureOnly
allSpeciesSite <- invert_2_ReadsOnly


# Site x Species table
cocoaSiteSpecies <- invert_2_cocoaSiteSpecies
pastureSiteSpecies <- invert_2_pastureSiteSpecies
allSiteSpecies <- invert_2_siteSpecies


allSiteSpecies.landUse <- ifelse(grepl("[0-9][0-9]C",rownames(allSiteSpecies)),
                                 yes = "Cocoa",
                                 no = ifelse(grepl("[0-9][0-9]P", rownames(allSiteSpecies)),
                                             yes = "Pasture",
                                             no = "Forest"))



rm(list = ls()[grepl(ls(), pattern = "invert_2|_20")])


## ----- Create the compositional matrices ----------------------
# Just do this for the 2g matrices that will be used for further analysis.

zPatterns(allSiteSpecies[,], 0)
zPatterns(allSpeciesSite[,-1], 0)

# square-root Bayesian-multiplicative replacement of zeros with the cmultRepl()
# function (Ladin et al., 2021) 

# Most other code in GitHub uses CZM. e.g. See
# https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-Biplot
# and https://raw.githubusercontent.com/ggloor/CoDaSeq/6ff864aade46cd3c8b0eff3bb54d5460775f92cd/CoDaSeq/vignettes/CoDaSeq_vignette.Rnw
# This latter contends that this is the most principled method.

allSiteSpecies_0repl <- cmultRepl(allSiteSpecies, label = 0,
                                  method = "CZM")


#boxplot(aplus(X = allSiteSpecies_0repl[ , ], parts = colnames(allSiteSpecies[ , ])))

allSiteSpecies_comps <- cdt.acomp(x = allSiteSpecies_0repl) %>% 
    as_tibble(., rownames = NA)

## see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/

#boxplot(cdt.acomp(x = allSiteSpecies_0repl, ), dots = T)

remove(allSiteSpecies_0repl)



## ----- Biodiversity Indicator 4 -------------------------------
## Number of keystone/priority species due to intervention.

# setup
cocoaKSPresent <-
    cocoaKeySpecies %>%
    dplyr::select(which(colSums(.) > 0)) %>% colnames()
pastureKSPresent <-
    pastureKeySpecies %>%
    dplyr::select(which(colSums(.) > 0)) %>% colnames()
forestKSPresent <-
    forestKeySpecies %>%
    dplyr::select(which(colSums(.) > 0)) %>% colnames()
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

## Which species are in neither? e.g. due to filtering.
keyInNeither <- keyInvertList[!(keyInvertList %in% c(cocoaKSPresent, pastureKSPresent))]

## Which species are in forests Only
keyInForestOnly <- forestKSPresent[forestKSPresent %in% keyInNeither]

keyInNone2g <- keyInNeither[!(keyInNeither %in% forestKSPresent)]



## ----- Indicator species analysis --------------------------- 
## Indicator species analysis: which species are particularly associated with
## intervention sites versus counterfactual sites?

# Doing this without the forest ones.



temp <- allSiteSpecies_comps %>%
    dplyr::select(which(colnames(allSiteSpecies_comps) %in% keyInvertList))

indicatorSpecies <- repeat.multipatt(
    matrix.name = temp,
    cluster.name = allSiteSpecies.landUse,
    quiet = F,
    func = "indval.g",
    p.cutoff = 0.05
)

remove(temp)

# indicatorSpecies <-
#     indicatorSpecies[order(indicatorSpecies$groupname,
#                            indicatorSpecies$mean.stat,
#                            decreasing = TRUE),] %>%
#     dplyr::select(groupname,
#                   species,
#                   frequency.sp,
#                   mean.stat,
#                   min.p.val,
#                   max.p.val,
#     )





## ----- Question 1.8A ------------------------------------------
## Does the number of keystone/priority species (i.e. key pollinator species)
## change over time compared to counterfactual?

# Since we have a repeated measures design we need to use GLMM.

specnumber(cocoaKeySpecies)

# First create a number of keystone/priority species for each 

cocoaKeySpecies$updatedPlot <- rownames(cocoaKeySpecies)
cocoaKeySpecies <-
    cocoaKeySpecies %>% left_join(dplyr::select(siteLookup, c(updatedPlot, site, system)))

pastureKeySpecies$updatedPlot <- rownames(pastureKeySpecies)
pastureKeySpecies <-
    pastureKeySpecies %>% left_join(dplyr::select(siteLookup, c(updatedPlot, site, system)))


# We're using reads as abundance, which is closer to what the indicator is
# asking for than the compositional approach even if correlated and not equal to
# abundance.
tbSpeciesChange <- tibble(
    siteNames = c(cocoaKeySpecies$updatedPlot, pastureKeySpecies$updatedPlot),
    siteField = c(cocoaKeySpecies$site, pastureKeySpecies$site),
    siteType = c(cocoaKeySpecies$system, pastureKeySpecies$system),
    keySpeciesCount = c(specnumber(dplyr::select(
        cocoaKeySpecies, where(is.numeric)
    )),
    specnumber(dplyr::select(
        pastureKeySpecies, where(is.numeric)
    ))),
    keySpeciesAbundance = c(rowSums(dplyr::select(
        cocoaKeySpecies, where(is.numeric)
    )), rowSums(dplyr::select(
        pastureKeySpecies, where(is.numeric)
    )))
)


tbSpeciesChange$siteField <-
    factor(tbSpeciesChange$siteField,
           levels = unique(tbSpeciesChange$siteField))

temp1 <- tbSpeciesChange %>% group_by(siteField) %>% 
            summarise(meanCount = mean(keySpeciesCount))
temp2 <- tbSpeciesChange %>% group_by(siteField) %>% 
            summarise(meanAbundance = mean(keySpeciesAbundance))
temp2$x <- seq(0.75, nrow(temp2) + 0.25, by = 1)
temp2$xend <- seq(1.25, nrow(temp2) + 0.25, 1)


tbSpeciesChange <- left_join(tbSpeciesChange, temp1) %>% left_join(temp2)
tbSpeciesChange$siteField <- factor(tbSpeciesChange$siteField, levels = temp1$siteField)

remove(temp1, temp2)


library(lmerTest) # this may mask rand.

# keySpeciesCountLMER <- lme4::lmer(keySpeciesCount ~
#                                       siteType +
#                                       (1 | siteField),
#                                   data = tbSpeciesChange,
#                                   REML = TRUE)
# anova(keySpeciesCountLMER)
# lmerTest::rand(keySpeciesCountLMER)
# summary(keySpeciesCountLMER)
# test_keySpeciesCountLMER <- car::Anova(keySpeciesCountLMER)
# 
# 
# sjPlot::plot_model(keySpeciesCountLMER,
#                    show.values = T,
#                    show.p = T)


# Create graph where x-axis has each field, with individual points for each
# plot. y-axis should be number of species. Add in means for fields, means & 95%
# CI for cocoa/pasture, and possibly significance as well.

# 
# plotTBSpeciesCount <- tbSpeciesChange %>%
#     ggplot() +
#     geom_boxplot(aes(siteField, keySpeciesCount),
#                  outlier.shape = NA) +
#     geom_jitter(
#         aes(siteField, keySpeciesCount, color = siteType),
#         width = 0.1,
#         height = 0
#     ) +
#     # geom_segment(aes(
#     #     x = x,
#     #     y = meanCount,
#     #     xend = xend,
#     #     yend = meanCount
#     # ),
#     # color = "black",
#     # size = 1.25) +
#     #
#     # facet_grid(cols=vars(siteType), space = "free", scales = "free_x") +
#     labs(color = "Field Type", y = "Key Species Richness", x = "Field Site Code")




## ----- BI5 / Question 1.8B -------------------------------------------
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
    geom_boxplot(aes(siteType, keySpeciesAbundance),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, keySpeciesAbundance, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(y = "Key Species Abundance") +
    scale_x_discrete(name = "Field Type",
                     labels=c("COCOA" = "Shaded Cocoa",
                              "PASTURE" = "Pasture"
                              )) +
    scale_y_log10(oob = scales::squish_infinite) +
    theme(legend.position="")


#https://cmdlinetips.com/2019/05/how-to-highlight-select-data-points-with-ggplot2-in-r/
#https://ggplot2.tidyverse.org/reference/geom_segment.html

## _____ General species richness, abundance, & diversity questions ____________________

## ----- BI6 / Question 1.5 --------------------------------------------

## Does species richness change over time compared to counterfactual? &
## Biodiversity Indicator 6: Change in species richness due to interventions.
## H0: Increase in species richness in cocoa fields compared with pasture. H1:
## No change or decrease in species richness.

# note this differs from the above because it is all species richness, not just
# key species richness.

# Calculate species richness for cocoa and pasture, then add it to the existing
# tbSpeciesChange table

# Calculate Species Richness for cocoa
cocoaSpeciesRichness <- specnumber(cocoaSiteSpecies , MARGIN = 1)

mean(cocoaSpeciesRichness)
sd(cocoaSpeciesRichness)

# Calculate Species Richness for pasture
pastureSpeciesRichness <- specnumber(pastureSiteSpecies , MARGIN = 1)

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
    col = alpha(rgb(1,0,0), 0.5),
    lwd = 2,
    ci.lty = 0,
    ci.col = alpha(rgb(1,0,0), 0.25),
    main = "Species Accumulation",
    ylab = "Expected (mean) species richness"
)

plot(
    cocoaSpecAccum,
    ci.type = "poly",
    col = alpha(rgb(0,0,1), 0.5),
    lwd = 2,
    ci.lty = 0,
    ci.col = alpha(rgb(0,0,1), 0.25),
    add = TRUE
)

# now plot by field


tbSpeciesAllCount <- tbSpeciesChange %>%
    ggplot() +
    geom_boxplot(aes(siteType, allSpeciesRichness),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, allSpeciesRichness, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(y = "Species Richness") +
    scale_x_discrete(name = "Field Type",
                     labels=c("COCOA" = "Shaded Cocoa",
                              "PASTURE" = "Pasture"
                     )) +
    scale_y_log10(oob = scales::squish_infinite) +
    theme(legend.position="")


## ----- Question 1.6 ------------------------------
# NP Test for Question 1.6: Does relative abundance of species change over time
# compared to counterfactual? H0: Increase in relative abundance of species on
# cocoa fields compared with pastures. H1: No change or decrease in relative
# abundance of species.

# Calculate species richness for cocoa and pasture, then add it to the existing
# tbSpeciesChange table

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
    geom_boxplot(aes(siteType, allAbundance),
                 outlier.shape = NA) +
    geom_jitter(aes(siteType, allAbundance, color = siteType),
                width = 0.1,
                height = 0) +
    labs(color = "Field Type", y = "All Species Abundance (reads)")+
scale_x_discrete(name = "Field Type",
                 labels=c("COCOA" = "Shaded Cocoa",
                          "PASTURE" = "Pasture"
                 )) +
    scale_y_log10(oob = scales::squish_infinite) +
    theme(legend.position="")




## ----- BI7 / Question 1.7 --------------------------------------------
# Biodiversity Indicator 7: Change in biodiversity indices
# due to interventions. 
# 1.7:Are there changes in Shannon’s diversity index (and others) over time compared
# to counterfactual?

# Here, we'll want to calculate the site-specific diversity metrics. Hill
# numbers include species richness (Hill q = 0), effective species richness
# (Hill q = 1), and inverse Simpson index (Hill q = 2). Can cite Jost (2006).

diversityMetrics <- tibble(
    siteNames = rownames(allSiteSpecies),
    spRichness = vegan::specnumber(allSiteSpecies),
    shannon = vegan::diversity(allSiteSpecies,
                               index = "shannon",
                               MARGIN = 1),
    simpson =  vegan::diversity(allSiteSpecies,
                                index = "simpson",
                                MARGIN = 1),
    chao1 =  vegan::estimateR(allSiteSpecies)[2, ],
    effectiveSR = exp(shannon),
    invSimpson = vegan::diversity(allSiteSpecies,
                                  index = "invsimpson",
                                  MARGIN = 1)
)

tbSpeciesChange <- left_join(tbSpeciesChange, diversityMetrics)


# Test to compare site diversities between land use types
# Test for shannons.

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
    geom_boxplot(aes(siteType, shannon),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, shannon, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Shannon Diversity") +
    scale_x_discrete(name = "Field Type",
                     labels=c("COCOA" = "Shaded Cocoa",
                              "PASTURE" = "Pasture"
                     )) +
    #scale_y_log10(oob = scales::squish_infinite) +
    theme(legend.position="")


# test for Simpson's

simpsonLMER <- lme4::lmer(simpson ~
                              siteType +
                              (1 | siteField),
                          data = tbSpeciesChange,
                          REML = TRUE)
anova(simpsonLMER)
lmerTest::rand(simpsonLMER)
summary(simpsonLMER)
test_shannonLMER <- car::Anova(simpsonLMER)

tbSpeciesChange %>%
    ggplot() +
    geom_boxplot(aes(siteType, simpson),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, simpson, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Simpson Diversity") +
    scale_x_discrete(name = "Field Type",
                     labels=c("COCOA" = "Shaded Cocoa",
                              "PASTURE" = "Pasture"
                     )) +
    #scale_y_log10(oob = scales::squish_infinite) +
    theme(legend.position="")




## ----- Proposed Indicator 1: Alpha -------------------------------------------
# PI1: Alpha diversity

# Test to compare site diversities between land use types

# Test for ESR

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
    geom_boxplot(aes(siteType, effectiveSR),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, effectiveSR, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Effective Species Richness") + 
    scale_x_discrete(name = "Field Type",
                     labels=c("COCOA" = "Shaded Cocoa",
                              "PASTURE" = "Pasture"
                     )) +
    #scale_y_log10(oob = scales::squish_infinite) +
    theme(legend.position="")

# Test for inverse Simpson

invSimpsonLMER <- lme4::lmer(invSimpson ~
                              siteType +
                              (1 | siteField),
                          data = tbSpeciesChange,
                          REML = TRUE)
anova(invSimpsonLMER)
lmerTest::rand(invSimpsonLMER)
summary(invSimpsonLMER)
test_shannonLMER <- car::Anova(invSimpsonLMER)

tbSpeciesChange %>%
    ggplot() +
    geom_boxplot(aes(siteType, invSimpson),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteType, invSimpson, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Inv. Simpson Diversity") +
    scale_x_discrete(name = "Field Type",
                     labels=c("COCOA" = "Shaded Cocoa",
                              "PASTURE" = "Pasture"
                     )) +
    #scale_y_log10(oob = scales::squish_infinite) +
    theme(legend.position="")


## ----- PI2: Beta diversity w/ Aitchison distance -----------------------------

# helper function
get_lower_tri <- function(inpMatrix){
    inpMatrix[upper.tri(inpMatrix, diag = T)]<- NA
    return(inpMatrix)
}


aitchisonPlot <- vegdist(allSiteSpecies_comps, method = "eucl", diag = F)
# aitchison distance uses the euclidian distance of the compositional data that
# has been center log transformed; see
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7811025/

min(aitchisonPlot)
max(aitchisonPlot)

#aggregate by site--since we can treat transformed data as euclidian, use mean 

aitchisonSite <- allSiteSpecies_comps %>%
    tibble::rownames_to_column(var = "updatedPlot") %>%
    left_join(siteLookup) %>%
    arrange(., siteEasy) %>%
    group_by(siteEasy) %>%
    summarise(across(GSFX_000000168:ZSFX_000026178, mean)) %>%
    tibble::column_to_rownames("siteEasy") %>%
    vegdist("euclid", diag = F, upper = F)

min(aitchisonSite)
max(aitchisonSite)



## create plots
## good guide here: https://datavizpyr.com/heatmap-from-matrix-using-ggplot2-in-r/
# and here http://sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

aitchisonSite %>%
    as.matrix() %>%
    get_lower_tri() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("plot1") %>%
    pivot_longer(-c(plot1),
                 names_to = "plot2",
                 values_to = "distance",
                 values_drop_na = T) %>%
    ggplot(aes(x = plot1, y = plot2, fill = distance)) + 
    geom_raster() +
    scale_fill_gradient(low = "blue", high = "orange",  
                         name="Aitchison\nDistance") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1)) +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        #axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.5, 0.7),
        legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))






## ----- PI3: Change in beta diversity -----------------------------------------

# Looking at the overall distance between the treatments. 


aitchisonTreatment <- allSiteSpecies_comps %>%
    tibble::rownames_to_column(var = "updatedPlot") %>%
    left_join(siteLookup) %>%
    group_by(system) %>%
    summarise(across(GSFX_000000168:ZSFX_000026178, mean)) %>%
    tibble::column_to_rownames("system") %>%
    vegdist("euclid", diag = F)

# Plot of Aitchison distance BETWEEN and WITHIN groups is what was plotted in
# Ladin et al. So cocoa-cocoa, cocoa-pasture, pasture-pasture.

aitchisonPlot %>%
    as.matrix() %>%
    get_lower_tri() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("plot1") %>%
    pivot_longer(-c(plot1),
                 names_to = "plot2",
                 values_to = "distance",
                 values_drop_na = T) %>%
    mutate(type = if_else(grepl("P", plot1) & grepl("P", plot2),
                         "Pasture-Pasture",
                         if_else(grepl("C", plot1) & grepl("C",plot2),
                                 "Shd. Cocoa-Shd. Cocoa",
                                 if_else(grepl("01F", plot1) | grepl("01F", plot2), 
                                 "Forest",
                                 "Shd. Cocoa-Pasture"))
                         )
           ) %>%
    dplyr::filter(type != "Forest") %>%
    ggplot(aes(y = distance, x = type)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.25) +
    geom_boxplot(aes(fill = type)) +
    theme(legend.position = "none") +
    labs(x = element_blank(),
         y = "Aitchison Distance (Plots)"
          )
    






## ----- PI 4: qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# create a column for pretty names.
invertData$pretty <- paste0(invertData$species, "_", str_sub(invertData$id, -6,-1))

# Create PCA for all plots
temp <- allSiteSpecies_comps %>%
    tibble::rownames_to_column(var = "updatedPlot") %>%
    left_join(siteLookup)

# rename the species columns to plot nicely

colnames(temp) <- plyr::mapvalues(colnames(temp),
                                  invertData$id,
                                  invertData$pretty)
    # ignore the "not present in x' message.


# create plot pcas
pca_plots <- temp %>%
    dplyr::select(Class_Arachnida_000168:Class_Arachnida_026178) %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaPlots <- fviz_pca_ind(
    pca_plots,
    geom.ind = "point",
    col.ind = temp$system,
    addEllipses = T,
    ellipse.type = "convex",
    legend.title = "Group"
)
# ggpubr::ggpar(viz_pcaPlots,
#               title = "Plots - PCA")

viz_pcaPlots_contrib <- fviz_contrib(pca_plots, choice = "ind", axes = 1:2)

fviz_pca_biplot(pca_plots,
                # Sites
                col.ind = temp$system,
                addEllipses = T,
                ellipse.type = "convex",
                label = "var",
                repel = T,
                max.overlaps = 5,
                alpha.var ="contrib")


# 
temp <- allSiteSpecies_comps %>%
    tibble::rownames_to_column(var = "updatedPlot") %>%
    left_join(siteLookup) %>%
    group_by(siteEasy) %>%
    summarise(across(GSFX_000000168:ZSFX_000026178, mean)) %>%
    tibble::column_to_rownames("siteEasy")
    

temp$system <- ifelse(grepl("Cocoa",rownames(temp)),
                      yes = "Shd. Cocoa",
                      no = ifelse(grepl("Pasture", rownames(temp)),
                                  yes = "Pasture",
                                  no = "Forest"))

pca_sites <- temp %>%
    dplyr::select(GSFX_000000168:ZSFX_000026178) %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaSites <- fviz_pca_ind(
    pca_sites,
    geom.ind = "point",
    col.ind = temp$system,
    addEllipses = T,
    ellipse.type = "convex",
    legend.title = "Group",
    repel = TRUE
    )

viz_pcaSites_contrib <- fviz_contrib(pca_sites, choice = "ind", axes = 1:2)

ggpubr::ggpar(viz_pcaSites,
              title = "Sites - PCA")





## ----- END OF CODE ---------------------------