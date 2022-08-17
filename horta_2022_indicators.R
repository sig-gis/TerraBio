# This analyis is based on soil samples collected in June at Horta and processed
# by EcoMol.

# This file take the eDNA species table provided by EcoMol and determines the
# TerraBio biodiversity indicators. There are also some species accumulation
# curves to verify sampling completeness.

## ----- Data ingestion ---------------------------------------

# libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(stringr)
library(compositions)
library(zCompositions)
library(iNEXT)

source("functions.R")
source("../../../RCode/R_Scripts/triplet_fixer.R")
source("horta_2022_data_processing.R")

rare = 50

## ----- Accumulation curves ----------------------------------


