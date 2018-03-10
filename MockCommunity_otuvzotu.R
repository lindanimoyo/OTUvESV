#March 5, 2018
#check out the mock communities from T1, T2, T3 runs in OTU v ZOTU


#Reset R's Brain
rm(list=ls())

#install.packages("plyr")
library(plyr )
###Community analyses
#install.packages("vegan")
library(vegan)
library("stringr")
library("plyr")

#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")

#upload all mock communities
mock_T1_otu <- read.csv("data/16S_usearch10/T1_otu_mockandnegs.csv")
mock_T1_zotu <- read.csv("data/16S_usearch10/T1_zotu_mockandnegs.csv")

mock_T2_otu <- read.csv("data/16S_usearch10/T2_otu_mockandnegs.csv")
mock_T2_zotu <- read.csv("data/16S_usearch10/T2_zotu_mockandnegs.csv")

mock_T3_otu <- read.csv("data/16S_usearch10/T3_otu_mockandnegs.csv")
mock_T3_zotu <- read.csv("data/16S_usearch10/T3_zotu_mockandnegs.csv")

#subtract out one of teh neg controls so function works
mock_T3_otu <- mock_T3_otu[ ,-3]
mock_T3_zotu <- mock_T3_zotu[ ,-3]

#make function to clean up
cleanupfunction <- function(mockdataframe){
  #remove negative control and consensus lineage columns
  mockdataframe2 <- mockdataframe[ ,-c(3,4)]
  #make OTU ids the row names
  row.names(mockdataframe2) <- mockdataframe2$X
  #find out the class
  class(mockdataframe2)
  #make it a data frame and remove the X column
  mockdataframe2 <- as.data.frame(mockdataframe2[ ,-1])
  #make the row names the X column gain
  row.names(mockdataframe2) <-mockdataframe$X
  #give the dataframe column the name mock
  names(mockdataframe2) <- "mock"
  #rarefy to 10,000 sequences 
  mockdataframe2_rare <- rrarefy(mockdataframe2,10000)
  #get rowsums
  rowSums(mockdataframe2_rare)
  #transform dataframe
  mockdataframe2_rare <- t(mockdataframe2_rare)
  #convert from matrix to dataframe
  mockdataframe2_rare <- as.data.frame(mockdataframe2_rare )
  #add consensus lineage now
  mockdataframe2_rare$Consensus.Lineage <-  mockdataframe$Consensus.Lineage
  #find zeroes from rarefaction
  zeroes <- which(mockdataframe2_rare$mock==0)
  #remove zeroes from OTU table
  mockdataframe2_rare_nozeroes<- mockdataframe2_rare[-zeroes, ]
  #look at dataframe
  mockdataframe2_rare_nozeroes
  #sort data frame from highest to lowest mock
  mockdataframe2_rare_nozeroes_sorted <- mockdataframe2_rare_nozeroes[order(-mockdataframe2_rare_nozeroes[1]),]
  return(mockdataframe2_rare_nozeroes_sorted)
}

#make T1 OTU
mock_T1_otu_cleaned <- cleanupfunction(mock_T1_otu)
nrow(mock_T1_otu_cleaned) #18

mock_T1_zotu_cleaned <- cleanupfunction(mock_T1_zotu)
nrow(mock_T1_zotu_cleaned) #26

#make T2 otu
mock_T2_otu_cleaned <- cleanupfunction(mock_T2_otu)
nrow(mock_T2_otu_cleaned) #22

mock_T2_zotu_cleaned <- cleanupfunction(mock_T2_zotu)
nrow(mock_T2_zotu_cleaned) #33

#make T3 otu
mock_T3_otu_cleaned <- cleanupfunction(mock_T3_otu)
nrow(mock_T3_otu_cleaned) #28

mock_T3_zotu_cleaned <- cleanupfunction(mock_T3_zotu)
nrow(mock_T3_zotu_cleaned) #48

mock_T1_otu_cleaned
mock_T1_zotu_cleaned

mock_T2_otu_cleaned
mock_T2_zotu_cleaned

mock_T3_otu_cleaned
mock_T3_zotu_cleaned

mocktable <- rbind(mock_T1_otu_cleaned,mock_T1_zotu_cleaned,mock_T2_otu_cleaned,mock_T2_zotu_cleaned,mock_T3_otu_cleaned,mock_T3_zotu_cleaned )
mocktable$Timepoint <- c(rep("T1",nrow(mock_T1_otu_cleaned)),rep("T1",nrow(mock_T1_zotu_cleaned)),
                         rep("T2",nrow(mock_T2_otu_cleaned)),rep("T2",nrow(mock_T2_zotu_cleaned)),
                         rep("T3",nrow(mock_T3_otu_cleaned)),rep("T3",nrow(mock_T3_zotu_cleaned)))

mocktable$OTUvZOTU <- c(rep("OTU",nrow(mock_T1_otu_cleaned)),rep("ESV",nrow(mock_T1_zotu_cleaned)),
                         rep("OTU",nrow(mock_T2_otu_cleaned)),rep("ESV",nrow(mock_T2_zotu_cleaned)),
                         rep("OTU",nrow(mock_T3_otu_cleaned)),rep("ESV",nrow(mock_T3_zotu_cleaned)))

write.csv(mocktable, "Figures/otuvzotu_taxonomy/mockcommunities_otuvzotu.csv")

mocktable$ID <- paste(mocktable$Timepoint, mocktable$OTUvZOTU, sep="_")
mocktable$OTU_ID <-  paste(row.names(mocktable), mocktable$ID, sep="_")

#taxa in mock https://d2gsy6rsbfrvyb.cloudfront.net/media/amasty/amfile/attach/_D6300_ZymoBIOMICS_Microbial_Community_Standard_v1.1.3.pdf
#taxa in mock: 
"Pseudomonas aeruginosa"
"  Escherichia coli"
"Salmonella enterica"
"Lactobacillus fermentum"
"Enterococcus faecalis" 
"Staphylococcus aureus"
"Listeria monocytogenes "
  "Bacillus subtilis"
  
  #
  otu_taxa1 <- ldply(str_split(string = mocktable$Consensus.Lineage, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
  names(otu_taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  otu_taxa2 <- as.data.frame(lapply(otu_taxa1, gsub, pattern=" ", replacement=""))
  
  
  # 4.1. Make dataframe for phylum composition
  otu_abundance_taxa3 <- cbind(otu_taxa2,mocktable )  #convert the rownames to a proper column of the data.frame
  
 rownames(otu_abundance_taxa3) <- NULL
  
  otu_abundance_taxa3$Phylum <- as.character(otu_abundance_taxa3$Phylum)  #convert empty column to unidentified
  otu_abundance_taxa3$Phylum[is.na(otu_abundance_taxa3$Phylum)]<-"p__unidentified"
  
  otu_abundance_taxa3$Family <- as.character(otu_abundance_taxa3$Family)  #convert empty column to unidentified
  otu_abundance_taxa3$Family[is.na(otu_abundance_taxa3$Family)]<-"f__unidentified"
  otu_abundance_taxa3$Family[otu_abundance_taxa3$Family=="f__"]<-"f__unidentified"
  
  
  otu_abundance_taxa3$Kingdom_Phylum <- str_c(otu_abundance_taxa3$Phylum,otu_abundance_taxa3$Family, sep=";") #combine kingdom and phylun
  otu_abundance_taxa4 <- otu_abundance_taxa3[, !(colnames(otu_abundance_taxa3) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))]
  
  head(otu_abundance_taxa4)
  row.names(otu_abundance_taxa4) <- otu_abundance_taxa4$OTU_ID
  head(otu_abundance_taxa4)
  
  #remove consensus lineage, time point, otuvzotu, otu_ID column
  myvars <- names(otu_abundance_taxa4) %in% c("Consensus.Lineage", "Timepoint", "OTUvZOTU", "OTU_ID")
  otu_abundance_taxa5  <-otu_abundance_taxa4[!myvars]
  otu_abundance_taxa5
  
  # 4.2. Caluculate mean value for each phylum on each date
  library("reshape2")
  otu_abundance_taxa6 <- melt(otu_abundance_taxa5, id.vars= c("ID","Kingdom_Phylum"),variable.name="mock", value.name="Abundance")
  head( otu_abundance_taxa6 )
  class(otu_abundance_taxa6$Abundance)
  #get total number of sequences per sample - 10,0000
  otu_abundance_taxa6$Abundance <- (otu_abundance_taxa6$Abundance/10000)*100 # Check your sequence number and convert to relative abundance (%)
  head(otu_abundance_taxa6)
  
  phylum_abundance <- ddply(otu_abundance_taxa6, c("Kingdom_Phylum","ID"), summarise, sum = sum(Abundance, na.rm=TRUE))

  head(  phylum_abundance)
  

  #calculate mean abundance of each Phylum family combo by each time point by esv vs otu
  phylum_abundance2 <- ddply(phylum_abundance, c("Kingdom_Phylum","ID"), summarise, # Caluculate mean value on each timepoint
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
  phylum_over1 <- subset(phylum_abundance2, mean>0.08 & Kingdom_Phylum !="p__unidentified;f__unidentified") # Subset the phylum with >5%
  
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
  #34 families in the mock
  
  #ok so this script adds up the mean abundance of every phyla in Site D, then 100-that number, is the other phyla
  phylum_abundance_below1 <- ddply(phylum_abundance_over1, c("ID"), summarise, # Subset the phylum with <1% and classify into other phylum
                                   Kingdom_Phylum = "p__;f__other families",
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
  phylum_abundance3$Phylum <- ldply(str_split(string=phylum_abundance3$Kingdom_Phylum, pattern=";f__"), rbind)[,2] # Divide a column using "_".
  
  #order by most abundant phyla
  phyla_sorted <- phylum_abundance3[order(-phylum_abundance3[,3]), ]
  phyla_sorted 
  # 4.4. Reorder the phylum from the largest abundance to the lowest
  length(levels(phylum_abundance3$Phylum))
  

#re order the Sites - give site full names
phyla_sorted <- phylum_abundance3[order(phylum_abundance3[,2]), ]
phyla_sorted 
table(phyla_sorted$ID)

                            
# Graph making
library("ggplot2")
library("scales")
library("gridExtra")

fam_barplot1<-
  ggplot(phyla_sorted, aes(x=ID, y=mean, fill=Phylum)) + 
  geom_bar(stat="identity",color="black", width=0.4) +
  scale_fill_brewer(palette = "Paired")+    #When you use scale_fill_brewer(palette = "Paired"), the maximum is 12, or 13 ("other phylum" is not colored)
  labs(x=" ", y="Relative abundance (%)", fill="Family") + theme_bw()+
  theme(axis.text.x=element_text(size=12,angle=30, hjust=1),
        axis.text.y=element_text(size=12),
        legend.text = element_text(size=12)
  )


fam_barplot1


write.csv(phyla_sorted, "Figures/otuvzotu_taxonomy/Mock_bacteria_phylasorted_otu.csv")

ggsave(fam_barplot1,filename=paste("Figures/otuvzotu_taxonomy/fambarplot_mock_bacteria.pdf"), width=7, height=4.5)




