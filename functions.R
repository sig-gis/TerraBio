# This file contains functions for TerraBio eDNA analysis.

## ----- Data manip. funs. | transpose ---------------------------

speciesSite <- function(inputOTU, rare = 0) {
    #input OTU should only be data, no text/ID fields.
    # inputs Species as columns, outputs Species as rows.
     temp <- inputOTU[ , colSums(inputOTU) > rare]
     temp <- as.data.frame(inputOTU)
     temp.col.names <- colnames(temp)
     temp <- as.data.frame(t(temp))
     rownames(temp) <- temp.col.names
     
     return(temp)
 }

siteSpecies <- function(inputOTU, rowNames, rare = 0) {
    #input OTU should only be data, no text/descriptions.
    #inputs Species as rows, outputs Species as columns
    temp <- as.data.frame(inputOTU, row.names = rowNames)
    temp <- as.data.frame(t(temp))
    colnames(temp) <- rowNames
    temp <- temp[rowSums(temp) > rare,]
    
    return(temp)
    
}

## ----- Plotting funs. | Aitchison -------------------------------
# helper function
get_lower_tri <- function(inpMatrix){
    inpMatrix[upper.tri(inpMatrix, diag = T)]<- NA
    return(inpMatrix)
}

aitHeatmap <- function(input){
graph <-
    input %>%
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

return(graph)   
}
