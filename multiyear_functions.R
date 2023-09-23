# Functions for handling multi-year datasets. These functions assume that ALL
# YEARS have been compared with the same database using the same classifier to
# avoid inter-year differences in taxonomic assignment.

# Last modified by KDyson 09/15/2023

## ----- Data manip. funs. | compile ---------------------------

# Compiles multiple years of data into one table. Year column is created when
# merging. Tables should otherwise have the same columns prior to merging.
# Double check that all columns with the same name have the same data and data
# type.

# years should be a list of all years
# ... should be the dataframes for all years
# These MUST be in the same order otherwise you'll mess up your data.

# table1 <- data.frame('species' = 1:10)
# table2 <- data.frame('species' = 15:20)
# years <- c(2002,3003)
# multiyear_merge(years, table1,table2)

multiyear_merge <- function(years, ...){
    
    # rbind the different tables to one another
    output <- rbind(...)
    
    # count rows of each table input
    i = 1
    n <- ...length()
    length_list <- rep(0,n)

    for(i in 1:n) {
        length_list[[i]] <- nrow(...elt(i))
        print(...elt(i))
    }
    
    print(...length())
    print(length_list)

    # add the year column
    
    output$year <- rep(years, length_list)
    
    return(output)
    
}


## ----- Data manip. funs. | trim ---------------------------


