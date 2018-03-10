#this function gets the last column of the otu data frame, and makes a data frame of just the taxonomy ID

getdd <- function(dataframe){
    lastcolumn = ncol(dataframe)
    taxon = dataframe[ ,lastcolumn]
    id =row.names(dataframe)
    id
    dd <- data.frame(id,taxon)
    names <- gsub('.*; s__', '', dd$taxon)
    dd$Species <- names
    dd
}


#make function for transformed OTU tables
gettrans <- function(dataframe){
  lastcolumn = ncol(dataframe)
  otu.trans = t(dataframe[ ,-lastcolumn])
}