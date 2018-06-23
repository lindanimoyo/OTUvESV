#June 22, 2018
#get dataframe for fungi ESV top genera over 1% relative abundance 

#Reset R's Brain
rm(list=ls())

## Set working directory
#setwd("")
#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")

#load in desired libraries
library("ggplot2")
library("scales")
library("gridExtra")
library(plyr )
library(tidyverse)
library(stringr)
library(ggpubr)
library(vegan)

#source in functions
source('~/Dropbox/StatsandProgramming/source/gettaxondd.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/getrowsums.R', chdir = TRUE)

################################################################################
############################## 1. DATA PROCESSING ##############################
################################################################################

zotu_abundance_taxa <- read.csv("data/fungi_ITS2/zotutab_inoculum_taxonomy_fungi.csv", row.names=1)

head(zotu_abundance_taxa )
# 1.2. Extract otu-abundance table from OTU-abundance-taxon table
zotu_abundance <- zotu_abundance_taxa[, colnames(zotu_abundance_taxa) != "Consensus.Lineage"]
head(zotu_abundance)

#transform
zotu_abundance_trans <- t(zotu_abundance)
#rarefy
set.seed(10)
getrowsums(zotu_abundance_trans)

zotu_abundance<- rrarefy(zotu_abundance_trans,6660)

rowSums(zotu_abundance)
#rarefied
zotu_abundance <- t(zotu_abundance)
colSums(zotu_abundance)

#head(otu_taxa)
dd <- as.data.frame(zotu_abundance_taxa$Consensus.Lineage)
dd$OTU <- row.names(zotu_abundance_taxa)
names(dd) <- c("taxon","id")

# 1.3. Extract row of taxonomy from OTU-abundance-taxon table and modify the format
#otu_taxa <- otu_abundance_taxa["ConsensusLineage"]  # Extract taxon table

library("stringr")
library("plyr")
zotu_taxa1 <- ldply(str_split(string = dd$taxon, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(zotu_taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
zotu_taxa2 <- as.data.frame(lapply(zotu_taxa1, gsub, pattern=" ", replacement=""))
zotu_taxa3<- cbind(dd[,1:2 ],zotu_taxa2)

#check column names of OTU table and taxonomy match
row.names(zotu_abundance)==(zotu_taxa3$id)

zotu_abundance_taxa2 <- cbind(zotu_abundance, zotu_taxa3) #Combine otu ID and taxon table
head(zotu_abundance_taxa2)

##########################################################################################
########################## 4. Community composition at genus level ######################
##########################################################################################


# 4.1. Make dataframe for phylum composition
zotu_abundance_taxa3 <- cbind(OTU_ID = rownames(zotu_abundance_taxa2), zotu_abundance_taxa2)  #convert the rownames to a proper column of the data.frame
rownames(zotu_abundance_taxa3) <- NULL

zotu_abundance_taxa3$Phylum <- as.character(zotu_abundance_taxa3$Phylum)  #convert empty column to unidentified
zotu_abundance_taxa3$Phylum[is.na(zotu_abundance_taxa3$Phylum)]<-"p__unidentified"

zotu_abundance_taxa3$Family <- as.character(zotu_abundance_taxa3$Family)  #convert empty column to unidentified
zotu_abundance_taxa3$Family[is.na(zotu_abundance_taxa3$Family)]<-"f__unidentified"
zotu_abundance_taxa3$Family[zotu_abundance_taxa3$Family=="f__"]<-"f__unidentified"

zotu_abundance_taxa3$Genus <- as.character(zotu_abundance_taxa3$Genus)  #convert empty column to unidentified
zotu_abundance_taxa3$Genus[is.na(zotu_abundance_taxa3$Genus)]<-"g__unidentified"
zotu_abundance_taxa3$Genus[zotu_abundance_taxa3$Genus=="g__"]<-"g__unidentified"

zotu_abundance_taxa3$Kingdom_Phylum <- str_c(zotu_abundance_taxa3$Phylum,zotu_abundance_taxa3$Family, zotu_abundance_taxa3$Genus, sep=";") #combine kingdom and phylun
zotu_abundance_taxa4 <- zotu_abundance_taxa3[, !(colnames(zotu_abundance_taxa3) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))]

zotu_abundance_taxa4 <- zotu_abundance_taxa4[, -c(22,23)]
head(zotu_abundance_taxa4)


# 4.2. Caluculate mean value for each phylum on each date
library("reshape2")
zotu_abundance_taxa5 <- melt(zotu_abundance_taxa4, id.vars= c("OTU_ID","Kingdom_Phylum"), variable.name="Sample", value.name="Abundance")
zotu_abundance_taxa5 <- subset(zotu_abundance_taxa5, Abundance>0)
#get total number of sequences per sample
totalseqs <- colSums(zotu_abundance_taxa4[ ,2:21])
avgnumseqs <- sum(totalseqs)/20

zotu_abundance_taxa5$Abundance <- (zotu_abundance_taxa5$Abundance/avgnumseqs)*100 # Check your sequence number and convert to relative abundance (%)
head(zotu_abundance_taxa5)

phylum_abundance <- ddply(zotu_abundance_taxa5, c("Kingdom_Phylum","Sample"), summarise,
                          sum = sum(Abundance, na.rm=TRUE)
)
head(phylum_abundance)

library("stringr")
library("plyr")

phylum_abundance$Site <- str_sub(phylum_abundance$Sample,2,2) #get a column of site names
head(phylum_abundance)


phylum_abundance2 <- ddply(phylum_abundance, c("Kingdom_Phylum","Site"), summarise, # Caluculate mean value on each date
                           mean = mean(sum, na.rm=TRUE),
                           sd = sd(sum, na.rm=TRUE),
                           n = sum(!is.na(sum)),
                           se = sd/sqrt(n),
                           max = max(sum, na.rm=TRUE) 
)
head(phylum_abundance2)


# 4.3. Classify phylum with mean abundance < XX% (1% here) into "other phylum". 
# Twelve phylum and others are recommended to use palette="Paired".
# "maximum abundance" rather than "mean abundance" will be better when compositional change is large.
phylum_over1 <- subset(phylum_abundance2, mean>1) # Subset the genus with >1%

#if one site has more than 1% of that phylum then it stays
phylum_over1_vector <- as.vector(unique(phylum_over1$Kingdom_Phylum)) # Subset the phylum with >1%
length(phylum_over1_vector) #37

write.csv(phylum_over1, "Figures/otuvzotu_taxonomy/fungi_esv_genusover1%relativeabundance.csv")
