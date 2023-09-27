require(showtext)

## Colorways from the Alliance Transition Branding GUidelines (August 2020) 

corporateColorPalette <- c("#2F8927", "#003580", "#CC3333", "#009933", "#414042", "#19191A")
supportingColorPalette <- c("#F68B33", "#F5D226", "#8EBF3F", "#CAD32B", "#0088C6", "#9D9FA2", "#6E4C1F", "#993399")


## Specify fonts, uncomment below for first run.

#install.packages('showtext', dependencies = TRUE)
library(showtext)
font_add_google("Open Sans")
font_add("Calibri", regular = "calibri.ttf", bold = "calibrib.ttf", italic = "calibrii.ttf", bolditalic = "calibriz.ttf")

# check fonts available
font_families()

#Calibri is the alliance's web approved font. Open Sans family is the print approved.
fontWeb <- "Calibri" 
fontPrint <- "Open Sans"
