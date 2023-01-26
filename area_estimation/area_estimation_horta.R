#####################################################################
library(survey)
library(tidyverse)
library(knitr)
library(rmarkdown)
library(tidyr)
library(networkD3)
library(tidyverse)
library(dplyr)
library(scales)
library(lubridate)

# setwd("C:\\Users\\karis\\Documents\\SBP\\AreaEstimation")
path <- "C:/Users/hayle/Desktop/terrabio/"

# preprocess data ####
## read in CEO data ####
dataCEOYR_TB <- read.csv(paste(path,'ceo-TerraBio_Validation_2017_2021_LandTrendrResults_Horta-sample-data-2022-11-29.csv', sep = "")) 
head(dataCEOYR_TB)
colnames(dataCEOYR_TB)

## read in and prep GEE data ####
dataGEE_TB <- read.csv(paste(path,'ceo-standage-horta-2017-2021-clean.csv', sep = ""))
colnames(dataGEE_TB)
colnames(dataGEE_TB)[14]<-'Loss17_21'  # "Forest.loss.2017.2021." T/F
colnames(dataGEE_TB)[13]<-'LC2021'  # T/F, T = forest, F = non-forest
colnames(dataGEE_TB)[12]<-'Gain17_21'
colnames(dataGEE_TB)[28] <- 'pl_strata' # pl_horta_results, changing to be consistent with old script
colnames(dataGEE_TB)[21] <- 'est_yr'
dataGEE_TB<-dataGEE_TB[,c(29, 28, 20, 22, 21, 14, 13, 12, 15, 6:10)]
colnames(dataGEE_TB)
head(dataGEE_TB)

## Prepare CEO data for merging with GEE data ####
colnames(dataCEOYR_TB)
colnames(dataCEOYR_TB)[16]<-'LC_forest_2021CEO'  # Yes/No
colnames(dataCEOYR_TB)[17]<-'Forest_Loss_17_21'  # Yes/No
colnames(dataCEOYR_TB)[21]<-'Forest_Gain_17_21' # Yes/No
colnames(dataCEOYR_TB)[18]<-'YrLossCEO'  # yyyy 9999 or NA
colnames(dataCEOYR_TB)[22]<-'YrGainCEO'  # yyyy 9999 or NA
colnames(dataCEOYR_TB)[24] <- 'notes'

dataCEOYR_TB$YrLossCEO[is.na(dataCEOYR_TB$YrLossCEO)] <- 9999
dataCEOYR_TB$YrGainCEO[is.na(dataCEOYR_TB$YrGainCEO)] <- 9999

## merge GEE & CEO data, remove extra strata ####
dataall_TB <- merge(dataGEE_TB, dataCEOYR_TB[,c(5, 12, 16, 18, 22)],
                by.x = c('pl_sampleid', 'email'), by.y = c('pl_sampleid', 'email'))
colnames(dataall_TB)  # "LC_forest_2021CEO""YrLossCEO""YrGainCEO" from dataCEOYR

## read in and prep strata pixel-count data ####
dataStrata_TB <- read.csv(paste(path, "countsReadable_Horta_lt_loss_greatest_2017_2021_SAMZguidance_Neutral_v2.csv", sep = ""))
head(dataStrata_TB)
dataStrata_TB[,c(3, 2, 4)]
dataStrata_TB[4]<-'strataName'  # readable column now filled with 'strataName'
# readable column used to contain Forest, Degradation, Deforestation, NonForest

## Merge strata pixel-count data with CEO-GEE-combo data ####
#CEO has pl_strata
dataSBP_TB<- merge(dataall_TB, dataStrata_TB[c(3, 2, 4)], by.x= 'pl_strata', by.y = 'map_value', all.x = T)
head(dataSBP_TB)
colnames(dataSBP_TB)  # count, readable from dataStrata
rm(dataCEOYR_TB, dataall_TB, dataGEE_TB, dataStrata_TB)

######################

## convert CEO info to stand age ####
# goal: fill this df w/ stand ages (a row for each plot)
CEOStandAge_TB <- data.frame(matrix(ncol = length(2021:2017),
                                 nrow = length(dataSBP_TB$LC_forest_2021CEO)))
colnames(CEOStandAge_TB) <- 2021:2017

#### preprocess CEO info in dataSBP ####
# add loss & gain event T/F indicators
dataSBP_TB$HasLoss <- dataSBP_TB$YrLossCEO != 9999  # has forest loss 17-21 or not
dataSBP_TB$HasGain <- dataSBP_TB$YrGainCEO != 9999  # has forest gain 17-21 or not

# change columns to T/F
LC_forest_2021CEO_TF <- dataSBP_TB$LC_forest_2021CEO == 'Yes'
dataSBP_TB$LC_forest_2021CEO <- LC_forest_2021CEO_TF
flagged_TF <- dataSBP_TB$flagged == 'true'
dataSBP_TB$flagged <- flagged_TF
Loss17_21_TF <- dataSBP_TB$Loss17_21 == 'true'
dataSBP_TB$Loss17_21 <- Loss17_21_TF
LC2021_TF <- dataSBP_TB$LC2021 == 'true'
dataSBP_TB$LC2021 <- LC2021_TF
Gain17_21_TF <- dataSBP_TB$Gain17_21 == 'true'
dataSBP_TB$Gain17_21 <- Gain17_21_TF

# add info about loss & gain order 
dataSBP_TB$LossThenGain <- FALSE
dataSBP_TB$GainThenLoss <- FALSE
gain_and_loss_TB <- which(dataSBP_TB$HasLoss == TRUE & dataSBP_TB$HasGain == TRUE)
for (x in gain_and_loss_TB){
  if (dataSBP_TB[x, 'YrLossCEO'] > dataSBP_TB[x, 'YrGainCEO']) {
    dataSBP_TB$GainThenLoss[x] <- TRUE
  } else {
    dataSBP_TB$LossThenGain[x] <- TRUE
  }
}

#### filling stand age into CEOStandAge ####
StartAge_TB <- 47 # moderate assumed age of 15 in 1985

# returns a vector of stand ages from 2021 to 2017
# #StartAge is the assumed stand age in 2017 if that cannot be determined
# through other logic
# #special treatment if loss then gain: assume gain happened 1 year after loss,
# because gain is generally detected late
# #special treatment if the earliest known event is a forest gain,
# assume the stand age at the year of gain is 4 to account for the fact that
# forest gain is generally detected late
calc_stand_age_21_17 <- function(CEOinfo, StartAge_TB) {
  is_forest_2021 <- CEOinfo$LC_forest_2021CEO # T/F rather than 0 or 100
  yr_loss <- CEOinfo$YrLossCEO
  yr_gain <- CEOinfo$YrGainCEO
  has_loss <- CEOinfo$HasLoss
  has_gain <- CEOinfo$HasGain
  loss_then_gain <- CEOinfo$LossThenGain
  gain_then_loss <- CEOinfo$GainThenLoss
  # if is non-forest in 2021
  if (!is_forest_2021) {
    if (!has_loss & !has_gain) { # if not forest in 2021 w no loss or gain
      return(rep(0, 5)) #! value changes to 5 bc 2017-2021
      
    } else if (has_loss & !has_gain) { # if not forest in 2021 w loss and no gain
      stand_age_17_loss <- StartAge_TB:(StartAge_TB+(yr_loss-2017)-1)
      stand_age_since_loss <- rep(0, 2021-yr_loss+1)
      stand_age_17_21 <- c(stand_age_17_loss, stand_age_since_loss)
      return(rev(stand_age_17_21))
      
    } else if (!has_loss & has_gain) { # if not forest in 2021 w no loss and gain
      print('not possible 1')
      # assume is_forest_2021
      # assume stand age at year of gain is 4 already
      stand_age_since_gain <- (1+3):(2021-yr_gain+1+3)
      stand_age_bc_gain <- c(rep(0, 1000000), 1:3)
      stand_age_bc_21 <- c(stand_age_bc_gain,
                           stand_age_since_gain)
      return(rev(stand_age_bc_21)[1:5])
      
    } else if (has_loss & has_gain) { # if not forest in 2021 w loss and gain
      if (loss_then_gain) {
        print('not possible 2')
        # assume is_forest_2021
        # assume yr_gain to be 1 year after yr_loss
        yr_gain <- yr_loss + 1
        stand_age_since_gain <- 1:(2021-yr_gain+1)
        stand_age_from_loss_to_gain <- rep(0, yr_gain-yr_loss)
        stand_age_17_loss <- StartAge_TB:(StartAge_TB+(yr_loss-2017)-1)
        stand_age_17_21 <- c(stand_age_17_loss,
                             stand_age_from_loss_to_gain,
                             stand_age_since_gain)
        return(rev(stand_age_17_21))
        
        
      } else if (gain_then_loss) {
        stand_age_since_loss <- rep(0, 2021-yr_loss+1)
        # assume stand age at year of gain is 4 already
        stand_age_from_gain_to_loss <- (1+3):(yr_loss-yr_gain+3)
        stand_age_before_gain <- c(rep(0, 100000), 1:3)
        stand_age_bc_21 <- c(stand_age_before_gain,
                             stand_age_from_gain_to_loss,
                             stand_age_since_loss)
        return(rev(stand_age_bc_21)[1:5])
      } else {print("shouldn't be here!")}
      
    } else {print("shouldn't be here!")}
    
    # if is forest in 2021
  } else if (is_forest_2021) {
    if (!has_loss & !has_gain) { # if forest in 2021 with no loss or gain
      # return(MaxAge:(MaxAge-5+1))  # 47 to 43
      return((StartAge_TB+5-1):StartAge_TB)  # 51 to 80(StartAge)
      
    } else if (has_loss & !has_gain) { # if forest in 2021 w loss but no gain
      print('not possible 3')
      # assume !is_forest_2021
      stand_age_17_loss <- StartAge_TB:(StartAge_TB+(yr_loss-2017)-1)
      stand_age_since_loss <- rep(0, 2021-yr_loss+1)
      stand_age_17_21 <- c(stand_age_17_loss, stand_age_since_loss)
      return(rev(stand_age_17_21))
      
    } else if (!has_loss & has_gain) { #if forest in 2021 w gain and no loss
      # assume stand age at year of gain is 4 already
      stand_age_since_gain <- (1+3):(2021-yr_gain+1+3)
      stand_age_bc_gain <- c(rep(0, 1000000), 1:3)
      stand_age_bc_21 <- c(stand_age_bc_gain,
                           stand_age_since_gain)
      return(rev(stand_age_bc_21)[1:5])
      
    } else if (has_loss & has_gain) { # if forest in 2021 w loss and gain
      if (gain_then_loss) {
        print('not possible 4')
        # assume !is_forest_2021
        stand_age_since_loss <- rep(0, 2021-yr_loss+1)
        # assume stand age at year of gain is 4 already
        stand_age_from_gain_to_loss <- (1+3):(yr_loss-yr_gain+3)
        stand_age_before_gain <- c(rep(0, 100000), 1:3)
        stand_age_bc_21 <- c(stand_age_before_gain,
                             stand_age_from_gain_to_loss,
                             stand_age_since_loss)
        return(rev(stand_age_bc_21)[1:32])
        
      } else if (loss_then_gain) {
        # assume yr_gain to be 1 year after yr_loss
        yr_gain <- yr_loss + 1
        stand_age_since_gain <- 1:(2021-yr_gain+1)
        stand_age_from_loss_to_gain <- rep(0, yr_gain-yr_loss)
        stand_age_17_loss <- StartAge_TB:(StartAge_TB+(yr_loss-2017)-1)
        stand_age_17_21 <- c(stand_age_17_loss,
                             stand_age_from_loss_to_gain,
                             stand_age_since_gain)
        return(rev(stand_age_17_21))
        
      } else {print("shouldn't be here!")}
      
    } else {print("shouldn't be here!")}
  } else {print("shouldn't be here!")}
  
}

for (r in 1:nrow(dataSBP_TB)) {
  # print(r)
  CEOinfo_TB <- dataSBP_TB[r, c('LC_forest_2021CEO','YrLossCEO','YrGainCEO',
                          'HasLoss','HasGain','LossThenGain','GainThenLoss')]
  CEOStandAge_TB[r, ] <- calc_stand_age_21_17(CEOinfo_TB, StartAge_TB)
}
view(CEOStandAge_TB)

# #### compare with stand age and carbon by year from the map ####
# Carbon assumes a forest start age of 80 and is for natural regeneration
age_carbon_df_TB <- read.csv(paste(path, 'ceo-samples-standage-carbon-horta.csv', sep = ""))
age_df_TB <- age_carbon_df_TB[, grepl('standAge', colnames(age_carbon_df_TB))]
# sampleid column not exported, add it here
age_carbon_df_TB$sampleid <- c(3,6,8,9,11,16,19,26,28,43,46,47,53,59,69,73,48,85,86,87,81,82,83,84,0,1
,2,4,5,12,13,14,15,17,18,20,21,22,24,25,27,29,30,31,32,34,35,36,38,39,40,41,42,44,45,49,50,51,52,54,55,
56,57,58,60,62,63,64,65,66,67,68,70,71,72,74,76,79,7,10,23,33,37,61,75,77,78)
age_df_TB <- cbind(age_carbon_df_TB$sampleid, age_df_TB)
carbon_df_TB <- age_carbon_df_TB[, grepl('carbon', colnames(age_carbon_df_TB))]
carbon_df_TB <- cbind(age_carbon_df_TB$sampleid, carbon_df_TB)
# view(age_df[order(age_df[,1]),
#             order(colnames(age_df), decreasing = T)])
# view(CEOStandAge[order(dataSBP$pl_sampleid), ])
carbon_df_ord_TB <- carbon_df_TB[order(carbon_df_TB[,1]),
                           order(colnames(carbon_df_TB), decreasing = T)]

##### scatter plot carbon for each sample pt: carbon based on map (age_carbon_df) ####
# vs carbon based on CEO (CEOcarbonNReg from below made with start age 47)
# for 2021 only
plot(CEOcarbonNReg[order(dataSBP$pl_sampleid), 'carbon2021NReg'],
     carbon_df_ord$X2021_carbon*0.09,
     pch=19, xlab='based on CEO', ylab='based on map',
     main='Estimate carbon in sample plots (tonne) in 2021')
lines(0:7, 0:7)

## estimate carbon per pixel ###

### tropical humid forest south america: #### 
b0_TB <- 2.479989225520437
b1_TB <- 0.12579472538348987
b2_TB <- 5.180809811711873

### CI upper #! NEED TO UPDATE FOR TROP HUM SA VALUES
b0up_TB <- 83.6653
b1up_TB <- 0.048341
b2up_TB <- 1.286441

### CI lower: #! NEED TO UPDATE FOR TROP HUM SA VALUES
b0low_TB <- 62.1984
b1low_TB <- 0.086366
b2low_TB <- 4.068786

##CEO carbon calcs
CEOcarbonTropHum <- data.frame(matrix(ncol = 5, nrow = length(dataSBP_TB$LC_forest_2021CEO)))
for (i in 1:5){
  CEOcarbonTropHum[,i] <- b0_TB * (1-exp(-b1_TB * CEOStandAge_TB[i]) )^b2_TB * 0.09
  colnames(CEOcarbonTropHum)[i] <- paste0('carbon', 2022-i, 'TropHum')
}

##CEO carbon calcs
CEOcarbonTropHumUp <- data.frame(matrix(ncol = 5, nrow = length(dataSBP_TB$LC_forest_2021CEO)))
for (i in 1:5){
  CEOcarbonTropHumUp[,i] <- b0up_TB * (1-exp(-b1up_TB * CEOStandAge_TB[i]) )^b2up_TB * 0.09
  colnames(CEOcarbonTropHumUp)[i] <- paste0('carbon', 2022-i, 'TropHumUp')
}

##CEO carbon calcs
CEOcarbonTropHumLow <- data.frame(matrix(ncol = 5, nrow = length(dataSBP_TB$LC_forest_2021CEO)))
for (i in 1:5){
  CEOcarbonTropHumLow[,i] <- b0low_TB * (1-exp(-b1low_TB * CEOStandAge_TB[i]) )^b2low_TB * 0.09
  colnames(CEOcarbonTropHumLow)[i] <- paste0('carbon', 2022-i, 'TropHumLow')
}

#CarbonCon <- cbind(dataSBP, CEOcarbonCON, CEOcarbonCONlow, CEOcarbonCONup)
CarbonTropHum <- cbind(dataSBP_TB, CEOcarbonTropHum, CEOcarbonTropHumLow, CEOcarbonTropHumUp)

# view(round(CEOcarbonNReg[order(dataSBP$pl_sampleid), ],2))

rm(dataSBP_TB,CEOStandAge_TB, CEOcarbonTropHum, CEOcarbonTropHumLow, CEOcarbonTropHumUp)




# area weighted estimates of total carbon ####
### select data of interest - per-pixel carbon estimate and ci ####
C_est_ci_df_TB <- CarbonTropHum
forest_type_TB <- 'TropHum'
forest_type_full_TB <- 'tropical_humid'

### survey design ####
strat_design_TB <- svydesign(id = ~1, strata = ~pl_strata, fpc = ~count,
                          data = C_est_ci_df_TB)
strat_design_TB

### total carbon & confidence intervals based on ####
### 1) per-pixel carbon estimate, 2) upper CI of per-pixel C est
### and 3) lower CI of per-pixel C est.
C_est_ci_ci_df_TB <- data.frame()

for (yyyy in 2017:2021) { #! UPDATE FOR CORRECT YEARS
  # total C and CI based on per-pixel C estimate
  C_est_formula_TB <- as.formula(paste0('~carbon',as.character(yyyy),forest_type_TB))
  svy_tot_TB <- svytotal(C_est_formula_TB, strat_design_TB)  # total C
  ci_TB <- confint(svy_tot_TB)  # CI of total C
  C_est_ci_yyyy_TB <- cbind(coef(svy_tot_TB), ci_TB)
  # total C and CI based on upperCI of per-pixel C estimate
  C_est_formula_up_TB <- as.formula(paste0('~carbon',as.character(yyyy),
                                        forest_type_TB,'up'))
  svy_tot_up_TB <- svytotal(C_est_formula_up_TB, strat_design_TB)  # total C
  ci_up_TB <- confint(svy_tot_up_TB)  # CI of total C
  C_est_ci_yyyy_up_TB <- cbind(coef(svy_tot_up_TB), ci_up_TB)
  # total C and CI based on lowerCI of per-pixel C estimate
  C_est_formula_low_TB <- as.formula(paste0('~carbon',as.character(yyyy),
                                         forest_type_TB,'low'))
  svy_tot_low_TB <- svytotal(C_est_formula_low_TB, strat_design_TB)  # total C
  ci_low_TB <- confint(svy_tot_low_TB)  # CI of total C
  C_est_ci_yyyy_low_TB <- cbind(coef(svy_tot_low_TB), ci_low_TB)
  # save all results
  C_est_ci_ci_df_TB <- rbind(C_est_ci_ci_df_TB,
                          cbind(C_est_ci_yyyy_TB, C_est_ci_yyyy_up_TB, C_est_ci_yyyy_low_TB))
}
C_est_ci_ci_df_TB
C_est_ci_ci_df_TB <- cbind(2017:2021, C_est_ci_ci_df_TB)
colnames(C_est_ci_ci_df_TB) <- c('year',
                              'totalC_estimate_from_growthFunc_ton',
                              'upp95CI_totalC_estimate_from_growthFunc_ton',
                              'low95CI_totalC_estimate_from_growthFunc_ton',
                              'totalC_estimate_from_uppGrowthFunc_ton',
                              'upp95CI_totalC_estimate_from_uppGrowthFunc_ton',
                              'low95CI_totalC_estimate_from_uppGrowthFunc_ton',
                              'totalC_estimate_from_lowGrowthFunc_ton',
                              'upp95CI_totalC_estimate_from_lowGrowthFunc_ton',
                              'low95CI_totalC_estimate_from_lowGrowthFunc_ton')
view(C_est_ci_ci_df_TB)

write.csv(C_est_ci_ci_df_TB, file = paste0(path, 'results/',
                                        'C_estimate_95ci_loMiUpGrowthFunc_',
                                        'start',
                                        as.character(StartAge_TB),
                                        '_',
                                        forest_type_full_TB,
                                        '.csv'))
