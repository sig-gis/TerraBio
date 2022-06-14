## 2. Data analysis::Biodiversity 

## This file takes the processed eDNA data and conducts the biodiversity data
## analysis. This file specifically addresses the vertebrates using the Riaz
## primers.

## Note: here key species are threatened. all species are native species
## (non-domestic)

## There were two protocols followed, one using 2g of soil and the other using
## 20g. We used the 2g for analysis as it returned more reads.

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
library(dplyr)

rm(list = ls()[grepl(ls(), pattern = "inve|gill|zea")])

## ----- Summary of eDNA data ------------------------------

## eDNA reads found; total number of species variants, etc. 
riazData$total_reads <- rowSums(riazData[, grepl(x = colnames(riazData), "SFX[0-9]")])
# riazData_2$total_reads <- rowSums(riazData_2[, grepl(x = colnames(riazData_2), "SFX[0-9]")])
# riazData_20$total_reads <- rowSums(riazData_20[, grepl(x = colnames(riazData_2), "SFX[0-9]")])

totalReadsRiaz_unfiltered <- sum(riazData$total_reads)
totalReadsRiaz_filtered <- sum(riazData$total_reads[riazData$total_reads >= minAbun])
totalReadsRiaz_2 <- sum(riazData_2$total_reads)
totalReadsRiaz_20 <- sum(riazData_20$total_reads)

totalMOTURiaz_unfiltered <- length(unique(riazData$id[riazData$total_reads>0]))
totalMOTURiaz_filtered <- length(unique(riazData$id[riazData$total_reads>=minAbun]))
totalMOTURiaz_2 <- length(unique(riazData_2$id))
totalMOTURiaz_20 <- length(unique(riazData_20$id))



totalUniqueTaxonRiaz_unfiltered <- unique(riazData$family) %>%
    grep("Class_|Order_", ., invert = TRUE) %>% length()
totalUniqueTaxonRiaz_filtered <-
    unique(riazData$family[riazData$total_reads >= minAbun]) %>%
    grep("Class_|Order_", ., invert = TRUE) %>% length()
totalUniqueTaxonRiaz_2 <-
    unique(riazData_2$family) %>% grep("Class_|Order_", ., invert = TRUE) %>% length()
totalUniqueTaxonRiaz_20 <-
    unique(riazData_20$family) %>% grep("Class_|Order_", ., invert = TRUE) %>% length()

totalReadsRiaz_2 / totalReadsRiaz_filtered


## _______ Key species questions and indicators ________________

## ----- Key Species present -----------------------------------
# What keystone species are present in the samples?

# Cocoa key species
cocoaKeyVert_2 <- riazData_2_cocoaSiteSpecies[, colnames(riazData_2_cocoaSiteSpecies) %in% keyVertList]
# cocoaKeyVert_20 <- riazData_20_cocoaSiteSpecies[, colnames(riazData_20_cocoaSiteSpecies) %in% keyVertList]

cocoaKeyVert_2_Table <- tibble(
    id = colnames(cocoaKeyVert_2),
    `Cocoa Mean Abundance` = cocoaKeyVert_2 %>%
        lapply(., mean, MARGIN = 2) %>% 
        unlist() %>% round(4),
    `Cocoa Standard Deviation` = cocoaKeyVert_2 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)

# cocoaKeyVert_20_Table <- tibble(
#     id = colnames(cocoaKeyVert_20),
#     `Cocoa Mean Abundance` = cocoaKeyVert_20 %>%
#         lapply(., mean, MARGIN = 2) %>% 
#         unlist() %>% round(4),
#     `Cocoa Standard Deviation` = cocoaKeyVert_20 %>%
#         lapply(., sd) %>% 
#         unlist() %>% round(2)
# )


# Pasture key species
pastureKeyVert_2 <- riazData_2_pastureSiteSpecies[, colnames(riazData_2_pastureSiteSpecies) %in% keyVertList]
# pastureKeyVert_20 <- riazData_20_pastureSiteSpecies[, colnames(riazData_20_pastureSiteSpecies) %in% keyVertList]

pastureKeyVert_2_Table <-  tibble(
    id = colnames(pastureKeyVert_2),
    `Pasture Mean Abundance` = pastureKeyVert_2 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4),
    `Pasture Standard Deviation` = pastureKeyVert_2 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)

# pastureKeyVert_20_Table <-  tibble(
#     id = colnames(pastureKeyVert_20),
#     `Pasture Mean Abundance` = pastureKeyVert_20 %>%
#         lapply(., mean, MARGIN = 2) %>%
#         unlist() %>% round(4),
#     `Pasture Standard Deviation` = pastureKeyVert_20 %>%
#         lapply(., sd) %>% 
#         unlist() %>% round(2)
# )

# Forests
forestKeyVert_2 <-
    riazData_2_siteSpecies[grepl(x = rownames(riazData_2_siteSpecies), pattern = "SFX004-01F"), colnames(riazData_2_siteSpecies) %in% keyVertList]

# forestKeyVert_20 <-
#     riazData_20_siteSpecies[grepl(x = rownames(riazData_20_siteSpecies), pattern = "SFX004-01F"), colnames(riazData_20_siteSpecies) %in% keyVertList]
# all are 0

forestKeyVert_2_Table <-  tibble(
    id = colnames(forestKeyVert_2),
    `Forest Mean Abundance` = forestKeyVert_2 %>%
        lapply(., mean, MARGIN = 2) %>%
        unlist() %>% round(4),
    `Forest Standard Deviation` = forestKeyVert_2 %>%
        lapply(., sd) %>% 
        unlist() %>% round(2)
)
# forestKeyVert_20_Table <-  tibble(
#     id = colnames(forestKeyVert_20),
#     `Forest Mean Abundance` = forestKeyVert_20 %>%
#         lapply(., mean, MARGIN = 2) %>%
#         unlist() %>% round(4),
#     `Forest Standard Deviation` = forestKeyVert_20 %>%
#         lapply(., sd) %>%
#         unlist() %>% round(2)
# )

# create a table to present the data


keyVert_2_Table <- tibble(id = keyVertList) %>%
    left_join(keyVertSpecies[, c("id", "species")], "id") %>%
    left_join(cocoaKeyVert_2_Table, "id") %>%
    left_join(pastureKeyVert_2_Table, "id") %>%
    left_join(forestKeyVert_2_Table, "id")


# keyVert_20_Table <- tibble(
#     id = keyVertList
# ) %>%
#     left_join(keyVertSpecies[ , c("id", "species")], "id") %>%
#     left_join(cocoaKeyVert_20_Table, "id") %>%
#     left_join(pastureKeyVert_20_Table, "id") %>%
#     left_join(forestKeyVert_20_Table, "id")
    

remove(cocoaKeyVert_2_Table, cocoaKeyVert_20_Table,
       pastureKeyVert_2_Table, pastureKeyVert_20_Table,
       forestKeyVert_2_Table, forestKeyVert_20_Table)

## ----- Given what we know, use 2g going forward --------------

# Key species tables:
cocoaKeySpecies <- cocoaKeyVert_2
pastureKeySpecies <- pastureKeyVert_2
forestKeySpecies <- forestKeyVert_2



remove(cocoaKeyVert_2, cocoaKeyVert_20,
       pastureKeyVert_2, pastureKeyVert_20,
       forestKeyVert_2, forestKeyVert_20)

# Species x Site tables:
cocoaSpeciesSite <-
    riazData_2_cocoaOnly[riazData_2_cocoaOnly$id %in% vertList,]
pastureSpeciesSite <-
    riazData_2_pastureOnly[riazData_2_pastureOnly$id %in% vertList,]
forestSpeciesSite <-
    riazData_2_forestOnly[riazData_2_forestOnly$id %in% vertList,]
allSpeciesSite <-
    riazData_2_ReadsOnly[riazData_2_ReadsOnly$id %in% vertList,]

# Site x Species table
cocoaSiteSpecies <-
    riazData_2_cocoaSiteSpecies[, colnames(riazData_2_cocoaSiteSpecies) %in% vertList]
pastureSiteSpecies <-
    riazData_2_pastureSiteSpecies[, colnames(riazData_2_pastureSiteSpecies) %in% vertList]
allSiteSpecies <-
    riazData_2_siteSpecies[, colnames(riazData_2_siteSpecies) %in% vertList]

allSiteSpecies <- allSiteSpecies[ rowSums(allSiteSpecies) > 0, ]



allSiteSpecies.landUse <- ifelse(grepl("[0-9][0-9]C",rownames(allSiteSpecies)),
                                 yes = "Cocoa",
                                 no = ifelse(grepl("[0-9][0-9]P", rownames(allSiteSpecies)),
                                             yes = "Pasture",
                                             no = "Forest"))


rm(list = ls()[grepl(ls(), pattern = "riaz")])



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
    as_tibble(., rownames = "rownames") %>%
    tibble::column_to_rownames("rownames")


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
allKSPresent <- allSiteSpecies[, colnames(allSiteSpecies) %in% keyVertList] %>%
    colnames()

## Which key species are in cocoa and not in pasture? 
keyInCocoaOnly <- cocoaKSPresent[!(cocoaKSPresent %in% pastureKSPresent)]

## Which key species are in pasture and not in cocoa?
keyInPastureOnly <- pastureKSPresent[!(pastureKSPresent %in% cocoaKSPresent)]

## Which species are in both?
keyInBoth <- cocoaKSPresent[cocoaKSPresent %in% pastureKSPresent]

## Which species are in neither?
keyInNeither <- keyVertList[!(keyVertList %in% c(cocoaKSPresent, pastureKSPresent))]

## Which species are in forests Only
keyInForestOnly <- forestKSPresent[forestKSPresent %in% keyInNeither]

keyInNone2g <- keyInNeither[!(keyInNeither %in% forestKSPresent)]

## ----- Indicator species analysis --------------------------- 
## Indicator species analysis: which species are particularly associated with
## intervention sites versus counterfactual sites?


temp <- allSiteSpecies %>%
    dplyr::select(which(colnames(allSiteSpecies) %in% keyVertList))

indicatorSpecies <- repeat.multipatt(
    matrix.name = temp,
    cluster.name = allSiteSpecies.landUse,
    quiet = F,
    p.cutoff = 0.2
)

remove(temp)

# indicatorSpecies <-
#     indicatorSpecies[order(indicatorSpecies$groupname,
#                            indicatorSpecies$mean.stat,
#                            decreasing = TRUE),] %>%
#     select(groupname,
#            species,
#            frequency.sp,
#            mean.stat,
#            min.p.val,
#            max.p.val,
#     )







## ----- Question 1.8A ------------------------------------------
## Does the number of keystone/priority species (i.e. key pollinator species)
## change over time compared to counterfactual?

# Since we have a repeated measures design we need to use GLMM.

#specnumber(cocoaKeySpecies)

# First create a number of keystone/priority species for each 

cocoaKeySpecies$updatedPlot <- rownames(cocoaKeySpecies)
cocoaKeySpecies <-
    cocoaKeySpecies %>% left_join(
        dplyr::select(siteLookup, c(updatedPlot, site, system)))

pastureKeySpecies$updatedPlot <- rownames(pastureKeySpecies)
pastureKeySpecies <-
    pastureKeySpecies %>% left_join(
        dplyr::select(siteLookup, c(updatedPlot, site, system)))



tbSpeciesChange <- tibble(
    siteNames = c(cocoaKeySpecies$updatedPlot, pastureKeySpecies$updatedPlot),
    siteField = c(cocoaKeySpecies$site, pastureKeySpecies$site),
    siteType = c(cocoaKeySpecies$system, pastureKeySpecies$system),
    keySpeciesCount = c(specnumber(dplyr::select(
        cocoaKeySpecies, where(is.numeric)
    )),
    specnumber(
        dplyr::select(pastureKeySpecies, where(is.numeric))
    )),
    keySpeciesAbundance = c(rowSums(dplyr::select(
        cocoaKeySpecies, where(is.numeric)
    )), rowSums(dplyr::select(
        pastureKeySpecies, where(is.numeric)
    )))
)


tbSpeciesChange$siteField <-
    factor(tbSpeciesChange$siteField,
           levels = unique(tbSpeciesChange$siteField))

temp1 <-
    tbSpeciesChange %>% 
    group_by(siteField) %>% 
    summarise(meanCount = mean(keySpeciesCount))
temp2 <-
    tbSpeciesChange %>% 
    group_by(siteField) %>% 
    summarise(meanAbundance = mean(keySpeciesAbundance))
temp2$x <- seq(0.75, nrow(temp2) + 0.25, by = 1)
temp2$xend <- seq(1.25, nrow(temp2) + 0.25, 1)


tbSpeciesChange <- left_join(tbSpeciesChange, temp1) %>% 
    left_join(temp2)
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
# lmerTest::rand(keySpeciesCountLMER)
# 
# 
# sjPlot::plot_model(keySpeciesCountLMER, show.values=T, show.p = T)


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
    geom_boxplot(aes(siteField, keySpeciesAbundance),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, keySpeciesAbundance, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Key Species Abundance") +
    scale_x_discrete(name = "Field Site Code",
                     labels=c("SFX204-07C" = "Cocoa 1",
                              "SFX217-03C" = "Cocoa 2",
                              "SFX225-08C" = "Cocoa 3",
                              "SFX188-05C" = "Cocoa 4",
                              "SFX006-02C" = "Cocoa 5",
                              "SFX237-04P" = "Pasture 1",
                              "SFX051-02P" = "Pasture 2",
                              "SFX026-01P" = "Pasture 3",
                              "SFX128-07P" = "Pasture 4",
                              "SFX184-03P" = "Pasture 5"
                     )) +
    theme(legend.position="bottom")

    


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
lmerTest::rand(allSpeciesRichnessLMER)
summary(allSpeciesRichnessLMER)
test_allSpeciesRichnessLMER <- car::Anova(allSpeciesRichnessLMER)


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
    geom_boxplot(aes(siteField, allSpeciesRichness),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, allSpeciesRichness, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "All Species Richness") +
    scale_x_discrete(name = "Field Site Code",
                     labels=c("SFX204-07C" = "Cocoa 1",
                              "SFX217-03C" = "Cocoa 2",
                              "SFX225-08C" = "Cocoa 3",
                              "SFX188-05C" = "Cocoa 4",
                              "SFX006-02C" = "Cocoa 5",
                              "SFX237-04P" = "Pasture 1",
                              "SFX051-02P" = "Pasture 2",
                              "SFX026-01P" = "Pasture 3",
                              "SFX128-07P" = "Pasture 4",
                              "SFX184-03P" = "Pasture 5"
                     )) +
    theme(legend.position="bottom")


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
test_allSpeciesAbundanceLMER <- car::Anova(allAbundanceLMER)


# for specific species instead of totals; please add columns with that species' data.

# 
# tbSpeciesAllAbundance <- tbSpeciesChange %>%
#     ggplot() +
#     geom_boxplot(aes(siteField, allAbundance),
#                  outlier.shape = NA) +
#     geom_jitter(aes(siteField, allAbundance, color = siteType),
#                 width = 0.1,
#                 height = 0) +
#     labs(color = "Field Type", y = "All Species Abundance (reads)", x = "Field Site Code")
# 
# 



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
    geom_boxplot(aes(siteField, shannon),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, shannon, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Shannon Diversity") +
    scale_x_discrete(name = "Field Site Code",
                     labels=c("SFX204-07C" = "Cocoa 1",
                              "SFX217-03C" = "Cocoa 2",
                              "SFX225-08C" = "Cocoa 3",
                              "SFX188-05C" = "Cocoa 4",
                              "SFX006-02C" = "Cocoa 5",
                              "SFX237-04P" = "Pasture 1",
                              "SFX051-02P" = "Pasture 2",
                              "SFX026-01P" = "Pasture 3",
                              "SFX128-07P" = "Pasture 4",
                              "SFX184-03P" = "Pasture 5"
                     )) +
    theme(legend.position="bottom")


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
    geom_boxplot(aes(siteField, simpson),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, simpson, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Simpson Diversity") +
    scale_x_discrete(name = "Field Site Code",
                     labels=c("SFX204-07C" = "Cocoa 1",
                              "SFX217-03C" = "Cocoa 2",
                              "SFX225-08C" = "Cocoa 3",
                              "SFX188-05C" = "Cocoa 4",
                              "SFX006-02C" = "Cocoa 5",
                              "SFX237-04P" = "Pasture 1",
                              "SFX051-02P" = "Pasture 2",
                              "SFX026-01P" = "Pasture 3",
                              "SFX128-07P" = "Pasture 4",
                              "SFX184-03P" = "Pasture 5"
                     )) +
    theme(legend.position="bottom")









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
    geom_boxplot(aes(siteField, effectiveSR),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, effectiveSR, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Effective Species Richness") + 
    scale_x_discrete(name = "Field Site Code",
                     labels=c("SFX204-07C" = "Cocoa 1",
                              "SFX217-03C" = "Cocoa 2",
                              "SFX225-08C" = "Cocoa 3",
                              "SFX188-05C" = "Cocoa 4",
                              "SFX006-02C" = "Cocoa 5",
                              "SFX237-04P" = "Pasture 1",
                              "SFX051-02P" = "Pasture 2",
                              "SFX026-01P" = "Pasture 3",
                              "SFX128-07P" = "Pasture 4",
                              "SFX184-03P" = "Pasture 5"
                     )) +
    theme(legend.position="bottom")

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
    geom_boxplot(aes(siteField, invSimpson),
                 outlier.shape = NA) +
    geom_jitter(
        aes(siteField, invSimpson, color = siteType),
        width = 0.1,
        height = 0
    ) +
    labs(color = "Field Type", y = "Inv. Simpson Diversity") +
    scale_x_discrete(name = "Field Site Code",
                     labels=c("SFX204-07C" = "Cocoa 1",
                              "SFX217-03C" = "Cocoa 2",
                              "SFX225-08C" = "Cocoa 3",
                              "SFX188-05C" = "Cocoa 4",
                              "SFX006-02C" = "Cocoa 5",
                              "SFX237-04P" = "Pasture 1",
                              "SFX051-02P" = "Pasture 2",
                              "SFX026-01P" = "Pasture 3",
                              "SFX128-07P" = "Pasture 4",
                              "SFX184-03P" = "Pasture 5"
                     )) +
    theme(legend.position="bottom")


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
    summarise(across(RSFX_000000353:RSFX_000002791, mean)) %>%
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
                                  "Cocoa-Cocoa",
                                  if_else(grepl("01F", plot1) | grepl("01F", plot2), 
                                          "Forest",
                                          "Cocoa-Pasture"))
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




## ----- PI3: Change in beta diversity -----------------------------------------

# Looking at the overall distance between the treatments. 


aitchisonTreatment <- allSiteSpecies_comps %>%
    tibble::rownames_to_column(var = "updatedPlot") %>%
    left_join(siteLookup) %>%
    group_by(system) %>%
    summarise(across(RSFX_000000353:RSFX_000002791, mean)) %>%
    tibble::column_to_rownames("system") %>%
    vegdist("euclid", diag = F)







## ----- PI 4: qualitative assessment ------------------------------------------
library("FactoMineR")
library("factoextra")
# Good tutorial here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# Create PCA for all plots
temp <- allSiteSpecies_comps %>%
    tibble::rownames_to_column(var = "updatedPlot") %>%
    left_join(siteLookup)

pca_plots <- temp %>%
    dplyr::select(RSFX_000000353:RSFX_000002791) %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaPlots <- fviz_pca_ind(
    pca_plots,
    geom.ind = "point",
    col.ind = temp$system,
    addEllipses = T,
    ellipse.type = "convex",
    legend.title = "Group"
)
ggpubr::ggpar(viz_pcaPlots,
              title = "Plots - PCA")


temp <- allSiteSpecies_comps %>%
    tibble::rownames_to_column(var = "updatedPlot") %>%
    left_join(siteLookup) %>%
    group_by(siteEasy) %>%
    summarise(across(RSFX_000000353:RSFX_000002791, mean)) %>%
    tibble::column_to_rownames("siteEasy")


temp$system <- ifelse(grepl("Cocoa",rownames(temp)),
                      yes = "Cocoa",
                      no = ifelse(grepl("Pasture", rownames(temp)),
                                  yes = "Pasture",
                                  no = "Forest"))

pca_sites <- temp %>%
    dplyr::select(RSFX_000000353:RSFX_000002791) %>%
    PCA(., scale.unit = F, graph = F)

viz_pcaSites <- fviz_pca_ind(
    pca_sites,
    geom.ind = "point",
    col.ind = temp$system,
    addEllipses = T,
    ellipse.type = "convex",
    legend.title = "Group",
    mean.point = T
)
ggpubr::ggpar(pca_sites,
              title = "Sites - PCA")


## ----- END OF CODE ---------------------------