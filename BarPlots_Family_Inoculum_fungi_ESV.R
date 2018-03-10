#March 5, 2018


#Reset R's Brain
rm(list=ls())

## Set working directory
#setwd("")
#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")

#source in functions
source('~/Dropbox/StatsandProgramming/source/gettaxondd.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/getrowsums.R', chdir = TRUE)
################################################################################
############################## 1. DATA PROCESSING ##############################
################################################################################

# 1.1. Input the file of rarefied OTU-abundance-taxon table
#read in dataframes
#otu_abundance <- read.csv("data/fungi_ITS2/OTU_inoc_T1.csv", row.names=1)
zotu_abundance_taxa <- read.csv("data/fungi_ITS2/zotutab_inoculum_taxonomy.csv", row.names=1)

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

# 1.4. Figure out which taxonomic level to sum by - what percentage isnt identified to phyla or family?
#figure out what percentage of inoculum OTUs are only identified to kingdom bacteria
percentidentified <- function(otuabund){
  #calculate percentage with no Phyla ID - NA or "g__" for bacteria or "g__unidentified" for fungi
  totalphylanoID <- length(which(is.na(otuabund$Phylum))) + length(which(otuabund$Phylum == "p__")) + length(which(otuabund$Phylum == "p__unidentified"))
  percentwithphyla <- (nrow(otuabund) - totalphylanoID )/nrow(otuabund)*100
  
  #calculate percentage with no Class ID - NA or "g__" for bacteria or "g__unidentified" for fungi
  totalclassnoID <- length(which(is.na(otuabund$Class))) + length(which(otuabund$Class == "c__")) + length(which(otuabund$Class == "c__unidentified"))
  #calculation based on total number of rows minus rows with no ids/ total number of rows * 100
  percentwithclass <- (nrow(otuabund) - totalclassnoID )/nrow(otuabund)*100
  
  #calculate percentage with no Order ID - NA or "o__" for bacteria or "o__unidentified" for fungi
  totalordernoID <- length(which(is.na(otuabund$Order))) + length(which(otuabund$Order == "o__" )) + length(which(otuabund$Order == "o__unidentified"))
  #calculation based on total number of rows minus rows with no ids/ total number of rows * 100
  percentwithorder <- (nrow(otuabund) - totalordernoID )/nrow(otuabund)*100
  
  #calculate percentage with no Family ID - NA or "g__" for bacteria or "g__unidentified" for fungi
  totalfamnoID <- length(which(is.na(otuabund$Family))) + length(which(otuabund$Family == "f__")) + length(which(otuabund$Family == "f__unidentified"))
  #calculation based on total number of rows minus rows with no ids/ total number of rows * 100
  percentwithfamily <- (nrow(otuabund) - totalfamnoID)/nrow(otuabund)*100
  
  #calculate percentage with no Genus ID - NA or "g__" for bacteria or "g__unidentified" for fungi
  totalgenusnoID <- length(which(is.na(otuabund$Genus))) + length(which(otuabund$Genus == "g__")) + length(which(otuabund$Genus == "g__unidentified"))
  #calculation based on total number of rows minus rows with no ids/ total number of rows * 100
  percentwithgenus <- (nrow(otuabund) - totalgenusnoID)/nrow(otuabund)*100
  
  #calculate percentage with no species ID - NA or "s__" for bacteria or "s__unidentified" for fungi
  totalspeciesnoID <- length(which(is.na(otuabund$Species))) + length(which(otuabund$Species == "s__")) + length(which(otuabund$Species == "s__unidentified"))
  #calculation based on total number of rows minus rows with no ids/ total number of rows * 100
  percentwithspecies <- (nrow(otuabund) - totalspeciesnoID)/nrow(otuabund)*100
  
  #make names column
  names <- c("Phyla","Class","Order","Family","Genus","Species")
  #make percentages column 
  percentages <- c(percentwithphyla,percentwithclass, percentwithorder,percentwithfamily, percentwithgenus,percentwithspecies)
  as.data.frame(cbind(names,percentages))
}


#run function
fungipercentID_zotu <- percentidentified(zotu_abundance_taxa2)

fungipercentID_zotu


write.csv(fungipercentID_zotu, "Figures/otuvzotu_taxonomy/fungipercent_ID_zotu_rarefied.csv")
##########################################################################################
########################## 4. Community composition at phylum level ######################
##########################################################################################


# 4.1. Make dataframe for phylum composition
zotu_abundance_taxa3 <- cbind(OTU_ID = rownames(zotu_abundance_taxa2), zotu_abundance_taxa2)  #convert the rownames to a proper column of the data.frame
rownames(zotu_abundance_taxa3) <- NULL

zotu_abundance_taxa3$Phylum <- as.character(zotu_abundance_taxa3$Phylum)  #convert empty column to unidentified
zotu_abundance_taxa3$Phylum[is.na(zotu_abundance_taxa3$Phylum)]<-"p__unidentified"

zotu_abundance_taxa3$Family <- as.character(zotu_abundance_taxa3$Family)  #convert empty column to unidentified
zotu_abundance_taxa3$Family[is.na(zotu_abundance_taxa3$Family)]<-"f__unidentified"
#zotu_abundance_taxa3$Family[zotu_abundance_taxa3$Family=="f__"]<-"f__unidentified"

zotu_abundance_taxa3$Kingdom_Phylum <- str_c(zotu_abundance_taxa3$Kingdom, zotu_abundance_taxa3$Phylum,zotu_abundance_taxa3$Family, sep=";") #combine kingdom and phylun
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
phylum_over1 <- subset(phylum_abundance2, mean>6.1 & Kingdom_Phylum !="k__Fungi;p__unidentified;f__unidentified") # Subset the phylum with >5%

#if one site has more than 1% of that phylum then it stays
phylum_over1 <- as.vector(unique(phylum_over1$Kingdom_Phylum)) # Subset the phylum with >1%
length(phylum_over1)

#stick with over 6% bc its 12 families which is the max that can fit
phylum_abundance_over1 <- subset(phylum_abundance2, Kingdom_Phylum==phylum_over1[1]) # Subset the phylum with >1%
for(i in 2:length(phylum_over1)){ 
  phylum_abundance_over1 <-merge(phylum_abundance_over1, subset(phylum_abundance2, Kingdom_Phylum==phylum_over1[i]), all=T)
}
phylum_abundance_over1 

#figure out how many phyla are in the inoculum
famininoculum <- as.vector(unique(phylum_abundance2$Kingdom_Phylum))
length(famininoculum)
#117 families in the inoculum

#ok so this script adds up the mean abundance of every phyla in Site D, then 100-that number, is the other phyla
phylum_abundance_below1 <- ddply(phylum_abundance_over1, c("Site"), summarise, # Subset the phylum with <1% and classify into other phylum
                                 Kingdom_Phylum = "k__;p__;f__other families",
                                 mean = 100-sum(mean),
                                 sd = 0,
                                 n = mean(n),
                                 se= 0
)

#ok so what this says is that only ~1% of the reads in desert belong to other phyla
phylum_abundance_below1

 
phylum_abundance3 <- merge(phylum_abundance_over1, phylum_abundance_below1, all=T)
head(phylum_abundance3)


# 4.4. Reformat the date description
phylum_abundance3$Phylum <- ldply(str_split(string=phylum_abundance3$Kingdom_Phylum, pattern=";p__"), rbind)[,2] # Divide a column using "_".

#order by most abundant phyla
phyla_sorted <- phylum_abundance3[order(-phylum_abundance3[,3]), ]
phyla_sorted 
# 4.4. Reorder the phylum from the largest abundance to the lowest
length(levels(phylum_abundance3$Phylum))


#phylum_abundance3$Phylum <- factor(phylum_abundance3$Phylum, 
 #                                  levels=c("other families", "unidentified","Pleosporaceae",
  #                                "Pleosporales_fam_Incertae_sedis","Dothideomycetes_fam_Incertae_sedis",
   #                                "Dothideaceae","Herpotrichiellaceae","Mycosphaerellaceae","Dothioraceae",
     #                              "Rhytismataceae","Phaeosphaeriaceae","Rutstroemiaceae"))
                                   
#give shorter names
levels(phylum_abundance3$Phylum)[4] <- "Ascomycota;f__Dothideomycetes"
levels(phylum_abundance3$Phylum)[9] <- "Ascomycota;f__Pleosporales"

#re order the Sites - give site full names
phyla_sorted <- phylum_abundance3[order(phylum_abundance3[,2]), ]
phyla_sorted 
table(phyla_sorted$Site)

phyla_sorted$Sitefullnames <- c(rep("Desert",11),rep("Grassland",11),rep("Pine-Oak",11),rep("Subalpine",11),rep("Scrubland",11))

phyla_sorted
phyla_sorted$Sitefullnames <- factor(phyla_sorted$Sitefullnames, 
                                   levels = c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))                                             

# 4.5 Graph making
library("ggplot2")
library("scales")
library("gridExtra")

fam_barplot1<-
  ggplot(phyla_sorted, aes(x=Sitefullnames, y=mean, fill=Phylum)) + 
  geom_bar(stat="identity",color="black", width=0.4) +
  scale_fill_brewer(palette = "Paired")+    #When you use scale_fill_brewer(palette = "Paired"), the maximum is 12, or 13 ("other phylum" is not colored)
  labs(x="Site", y="Relative abundance (%) ESV", fill="Family") + theme_bw()+
  theme(axis.text.x=element_text(size=12,angle=30, hjust=1),
        axis.text.y=element_text(size=12),
        legend.text = element_text(size=12)
        )


fam_barplot1


phyla_sorted$otuvesv <- rep("ESV", nrow(phyla_sorted))
write.csv(phyla_sorted, "Figures/otuvzotu_taxonomy/Fungi_Inoculum_familysorted_zotu.csv")

ggsave(fam_barplot1,filename=paste("Figures/otuvzotu_taxonomy/fambarplot_inoculum_fungi_zotu.pdf"), width=7, height=4.5)


