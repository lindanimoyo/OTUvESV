#June 22, 2018
#make correlation figure for all FUNGAL genera with over 1% relative abundance

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
library(dplyr)

#################################################################
#############read in FUNGI esv vs OTU summarized by GENUS   
#################################################################

################################################################################
#load in data
################################################################################

#fungi top OTU genera
otu_fungi_genus <- read.csv("Figures/otuvzotu_taxonomy/fungi_OTU_genusover1%relativeabundance.csv", row.names=1)

#fungi top ESV genera
esv_fungi_genus <- read.csv("Figures/otuvzotu_taxonomy/fungi_esv_genusover1%relativeabundance.csv", row.names=1)

class(otu_fungi_genus$Kingdom_Phylum)
class(esv_fungi_genus$Kingdom_Phylum)

otu_fungi_genus$Kingdom_Phylum == esv_fungi_genus$Kingdom_Phylum

################################################################################
#join based on Kingdom_Phylum and Site columns
################################################################################

## First check the number of rows in these two files
nrow(otu_fungi_genus) #62
nrow(esv_fungi_genus) #65

names(esv_fungi_genus)
#rename columns for esv tableby otu v esv
names(esv_fungi_genus)[3]<-"mean_esv"
names(esv_fungi_genus)[4]<-"sd_esv"
names(esv_fungi_genus)[5]<-"n_esv"
names(esv_fungi_genus)[6]<-"se_esv"
names(esv_fungi_genus)[7]<-"max_esv"

#use dply to perform a full join: http://rpubs.com/williamsurles/293454
esvandotu_fungi_genus_full <- full_join(otu_fungi_genus, esv_fungi_genus, by = c("Kingdom_Phylum","Site"))
length(esvandotu_fungi_genus_full$Kingdom_Phylum)

esvandotu_fungi_genus_inner <- inner_join(otu_fungi_genus, esv_fungi_genus, by = c("Kingdom_Phylum","Site"))
length(esvandotu_fungi_genus_inner$Kingdom_Phylum)

################################################
#############plot figure
##############################################
sitecolors  <- c("red","orange","green","blue","purple")
sitenames <- c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine")
shapes <- c(17,3,16,15,18)

#re-order and re-label Site names
esvandotu_fungi_genus_full$Site <- factor(esvandotu_fungi_genus_full$Site,levels=c("D","W","G","P","S"), labels=c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))


p1 <- ggplot(esvandotu_fungi_genus_full, aes(x=mean, y=mean_esv, group=Site, col=Site, label=Kingdom_Phylum)) + geom_point(size=2) + 
  theme_bw()+
  geom_smooth(method = "lm", aes(col=Site)) + #add linear smoother
  labs(x="Genus Mean Relative Abundance OTU Fungi", y="Genus Mean Relative Abundance ESV Fungi", shape="Site") + #change y axis label 
  theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=12),  #make y axis tick sizes bigger
        axis.title=element_text(size=14))+ #make y axis label larger
  geom_errorbar(aes(ymin=mean_esv-se_esv, ymax=mean_esv+se_esv,width=.1))+#add in y standard error bars
  geom_errorbarh(aes(xmin=mean-se, xmax=mean+se)) +#add in x standard error bars
  scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual shapes and change the legend name
  #scale_shape_manual(values=shapes, labels=sitenames) + #add in manual shapes and change the legend name
  stat_cor(method="pearson",size=3.5) #add in pearson correlations

p1
p1 + geom_text(size =2)
p1 + facet_wrap(~Site)
p1 + facet_wrap(~Site)+ geom_text(size=1)
p1 + facet_wrap(~Kingdom_Phylum)

pdf("Figures/otuvzotu_taxonomy/Fungi_genus_correlations_bysite_allgeneraover1.pdf", width=8, height=8)
p1 + facet_wrap(~Site)
dev.off()

pdf("Figures/otuvzotu_taxonomy/Fungi_genus_correlations_allgeneraover1.pdf", width=8, height=8)
p1 
dev.off()

pdf("Figures/otuvzotu_taxonomy/Fungi_genus_bygenera_correlations_allgeneraover1.pdf", width=8, height=8)
p1 + facet_wrap(~Kingdom_Phylum)
dev.off()
