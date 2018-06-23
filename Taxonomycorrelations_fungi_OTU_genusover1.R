#June 22, 2018
#get dataframe for fungi OTU top genera over 1% relative abundance 

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

# 1.1. Input the file of rarefied OTU-abundance-taxon table
#read in dataframes
#otu_abundance <- read.csv("data/fungi_ITS2/OTU_inoc_T1.csv", row.names=1)
otu_abundance_taxa <- read.csv("data/fungi_ITS2/otutab_inoculum_taxonomy_fungi.csv", row.names=1)

head(otu_abundance_taxa )
# 1.2. Extract otu-abundance table from OTU-abundance-taxon table
otu_abundance <- otu_abundance_taxa[, colnames(otu_abundance_taxa) != "Consensus.Lineage"]
head(otu_abundance)

#transform
otu_abundance_trans <- t(otu_abundance)
#rarefy
set.seed(10)
getrowsums(otu_abundance_trans)

library(vegan)
otu_abundance<- rrarefy(otu_abundance_trans,6665)

rowSums(otu_abundance)
#rarefied
otu_abundance <- t(otu_abundance)
colSums(otu_abundance)

#head(otu_taxa)
dd <- as.data.frame(otu_abundance_taxa$Consensus.Lineage)
dd$OTU <- row.names(otu_abundance_taxa)
names(dd) <- c("taxon","id")

# 1.3. Extract row of taxonomy from OTU-abundance-taxon table and modify the format
#otu_taxa <- otu_abundance_taxa["ConsensusLineage"]  # Extract taxon table

library("stringr")
library("plyr")
otu_taxa1 <- ldply(str_split(string = dd$taxon, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(otu_taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu_taxa2 <- as.data.frame(lapply(otu_taxa1, gsub, pattern=" ", replacement=""))
otu_taxa3<- cbind(dd[,1:2 ],otu_taxa2)

#check column names of OTU table and taxonomy match
row.names(otu_abundance)==(otu_taxa3$id)

otu_abundance_taxa2 <- cbind(otu_abundance, otu_taxa3) #Combine otu ID and taxon table
head(otu_abundance_taxa2)

##########################################################################################
########################## 4. Community composition at phylum level ######################
##########################################################################################


# 4.1. Make dataframe for phylum composition
otu_abundance_taxa3 <- cbind(OTU_ID = rownames(otu_abundance_taxa2), otu_abundance_taxa2)  #convert the rownames to a proper column of the data.frame
rownames(otu_abundance_taxa3) <- NULL

otu_abundance_taxa3$Phylum <- as.character(otu_abundance_taxa3$Phylum)  #convert empty column to unidentified
otu_abundance_taxa3$Phylum[is.na(otu_abundance_taxa3$Phylum)]<-"p__unidentified"


otu_abundance_taxa3$Family <- as.character(otu_abundance_taxa3$Family)  #convert empty column to unidentified
otu_abundance_taxa3$Family[is.na(otu_abundance_taxa3$Family)]<-"f__unidentified"
otu_abundance_taxa3$Family[otu_abundance_taxa3$Family=="f__"]<-"f__unidentified"

otu_abundance_taxa3$Genus <- as.character(otu_abundance_taxa3$Genus)  #convert empty column to unidentified
otu_abundance_taxa3$Genus[is.na(otu_abundance_taxa3$Genus)]<-"g__unidentified"
otu_abundance_taxa3$Genus[otu_abundance_taxa3$Genus=="g__"]<-"g__unidentified"

otu_abundance_taxa3$Kingdom_Phylum <- str_c(otu_abundance_taxa3$Phylum,otu_abundance_taxa3$Family, otu_abundance_taxa3$Genus, sep=";") #combine kingdom and phylun
otu_abundance_taxa4 <- otu_abundance_taxa3[, !(colnames(otu_abundance_taxa3) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))]

otu_abundance_taxa4 <- otu_abundance_taxa4[, -c(22,23)]
head(otu_abundance_taxa4)

# 4.2. Caluculate mean value for each phylum on each date
library("reshape2")
otu_abundance_taxa5 <- melt(otu_abundance_taxa4, id.vars= c("OTU_ID","Kingdom_Phylum"), variable.name="Sample", value.name="Abundance")
otu_abundance_taxa5 <- subset(otu_abundance_taxa5, Abundance>0)
#get total number of sequences per sample
totalseqs <- colSums(otu_abundance_taxa4[ ,2:21])
avgnumseqs <- sum(totalseqs)/20

otu_abundance_taxa5$Abundance <- (otu_abundance_taxa5$Abundance/avgnumseqs)*100 # Check your sequence number and convert to relative abundance (%)
head(otu_abundance_taxa5)

phylum_abundance <- ddply(otu_abundance_taxa5, c("Kingdom_Phylum","Sample"), summarise,
                          sum = sum(Abundance, na.rm=TRUE)
)
head(phylum_abundance)


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
phylum_over1 <- subset(phylum_abundance2, mean>1) # Subset the phylum with >1%

#if one site has more than 1% of that phylum then it stays
phylum_over1_vector <- as.vector(unique(phylum_over1$Kingdom_Phylum)) # Subset the phylum with >1%
length(phylum_over1_vector) #34

write.csv(phylum_over1, "Figures/otuvzotu_taxonomy/fungi_OTU_genusover1%relativeabundance.csv")
