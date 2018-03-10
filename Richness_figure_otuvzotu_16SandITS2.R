#Dec 13, 2017

#make a combo figure with overall mean richness by site by time point for bacteria and fungi
#make one for just t3 and make one for all time points
#change ZOTU to ESV
#zOTUs

#Bacteria richness 
#use same codes as bac to make nice Figures and anovas
#Use alpha richness values from robert edgar usearch v10 commands, they are all rarefied to 10,000 sequences each
# Alpha diversity
#usearchv10 -alpha_div otutab_survey_e10000.txt -output alpha_otu_survey.txt
#usearchv10 -alpha_div zotutab_survey_e10000.txt -output alpha_zotu_survey.txt

#16S library
#Reset R's Brain
rm(list=ls())

#install.packages("plyr")
library(plyr )
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")


###############################################################################################################
######Load in data frames
###############################################################################################################
#read in fungi otu and zotu data frames with mean richness per site
fungi_otu <- read.csv("data/fungi_T1T2T3_meanrichnessbysite_otu.csv", row.names=1)
fungi_zotu <-  read.csv("data/fungi_T1T2T3_meanrichnessbysite_zotu.csv",  row.names=1)

#read in bacteria otu and zotu data frames with mean richness per site
bac_otu <- read.csv("data/bac_T1T2T3_site_all_otu.csv", row.names=1)
bac_zotu <- read.csv("data/bac_T1T2T3_meanrichnessbysite_zotu.csv",row.names=1)

#make factors correct order
fungi_otu$Sitenames <- factor(fungi_otu$Sitenames ,levels=c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))
fungi_zotu$Sitenames <- factor(fungi_zotu$Sitenames ,levels=c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))
bac_otu$Sitenames <- factor(bac_otu$Sitenames ,levels=c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))
bac_zotu$Sitenames <- factor(bac_zotu$Sitenames ,levels=c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))

###############################################################################################################
######Make figure function
##############################################################################################################

#includes an x axis
my.gg <- function(df, yname,nudgenumber) {
  ggplot(data=df, aes(x=Sitenames, y=mean, col=Sitenames, label=Tukeylabels)) + geom_point(size=2) +
    labs(x=" ", y=yname, col="Site") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme_bw()+ #change to black and white theme
    theme(strip.text.x = element_text(size = 12, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
          axis.text.y=element_text(size=10),  #make y axis tick sizes bigger
          axis.title.y=element_text(size=14), #make y axis label larger
          panel.spacing.x = unit(0,"null"))+ #removing spacing between panels
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_color_manual(values=c("red", "orange", "green","blue","purple"))+#add in manual colors for points/lines
    scale_x_discrete(labels= c("D","Sc","G","P","S"))+ #this change the x axis tick labels to just hte first letter instead of the entire word
    geom_text(nudge_y=nudgenumber) + #this is for the Tukey labels, makes them black and nudges them up on y axis so they aren't directly on top of point
    guides(col=FALSE) #remove legend for color
    #scale_y_continuous(position = "right") #put y axis on right side
}

#removes x axis and ticks etc
my.gg.noxaxis <- function(df, yname,nudgenumber) {
  ggplot(data=df, aes(x=Sitenames, y=mean, col=Sitenames, label=Tukeylabels)) + geom_point(size=2) +
    labs(y=yname, col="Site") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme_bw()+ #change to black and white theme
    theme(strip.text.x = element_text(size = 11, colour = "black"), #make T1, T2, T3 labels bigger
          axis.title.x = element_blank(),  #remove x axis title
          axis.text.x = element_blank(), #remove x axist text
          axis.ticks.x = element_blank(), #remove x axist title
          axis.text.y=element_text(size=10),  #make y axis tick sizes bigger
          axis.title.y=element_text(size=12), #make y axis label larger
          panel.spacing.x = unit(0,"null"))+ #removing spacing between panels
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_color_manual(values=c("red", "orange", "green","blue","purple"))+#add in manual colors for points/lines
    geom_text(nudge_y=nudgenumber) + #this is for the Tukey labels, makes them black and nudges them up on y axis so they aren't directly on top of point
    guides(col=FALSE) #remove legend for color
    #scale_y_continuous(position = "right") #put y axis on right side
}




####change labels from T1 T2 T3 to 6, 12, 18 months
#make a labeller function a la https://stackoverflow.com/questions/3472980/ggplot-how-to-change-facet-labels
timepointnames <- c('T1'="6 mos",'T2'="12 mos",'T3'="18 mos")


p1 <- my.gg(bac_otu, "Bacterial OTU Richness",40) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p1a <- my.gg.noxaxis(bac_otu, "Bacterial OTU Richness",40) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))

p2 <- my.gg(bac_zotu, "Bacterial ESV Richness",55) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p2a <- my.gg.noxaxis(bac_zotu, "Bacterial ESV Richness",55) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))

p3 <- my.gg(fungi_otu, "Fungal OTU Richness",10) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p3a <- my.gg.noxaxis(fungi_otu, "Fungal OTU Richness",10) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))

p4 <- my.gg(fungi_zotu, "Fungal ESV Richness",15) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p4a <- my.gg.noxaxis(fungi_zotu, "Fungal ESV Richness",15) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))

#arrange with axes
ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom", align="hv")

#arrange without x axes
ggarrange(p1a,p2a,p3a,p4a,ncol = 2, nrow = 2,  common.legend = TRUE, legend="bottom", align="hv")



ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2,align="hv") %>%
  ggexport(filename = "Figures/zotuvotu_richness/Meanrichnessbysite_3timepoints_fungiandbac_zotuvotu.pdf")

ggarrange(p1a,p2a,p3a,p4a,ncol = 2, nrow = 2,  common.legend = TRUE, legend="bottom", align="hv") %>%
  ggexport(filename = "Figures/zotuvotu_richness/Meanrichnessbysite_3timepoints_fungiandbac_zotuvotu_NOAXES.pdf")


#arrange just bacteria 
bacteria_otuvesv <- ggarrange(p1,p2, ncol = 2, common.legend = TRUE, legend="bottom", labels=c("A","B"))
bacteria_otuvesv

#arrange just fungi
fungi_otuvesv <- ggarrange(p3,p4, ncol = 2, common.legend = TRUE, legend="bottom", labels=c("A","B")) #add titles A and B to plots
fungi_otuvesv


pdf("Figures/zotuvotu_richness/Meanrichnessbysite_3timepoints_fungi_zotuvotu.pdf", height=3, width=6)
fungi_otuvesv
dev.off()

pdf("Figures/zotuvotu_richness/Meanrichnessbysite_3timepoints_bac_zotuvotu.pdf", height=3, width=6)
bacteria_otuvesv
dev.off()

###############################################################################################################
#######make figures based on elevation order
##############################################################################################################
#change order based on elevation
fungi_otu_elev <- fungi_otu
fungi_otu_elev$Sitenames <- factor(fungi_otu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))
fungi_zotu_elev <- fungi_zotu
fungi_zotu$Sitenames <- factor(fungi_zotu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))

bac_otu_elev <- bac_otu
bac_otu_elev$Sitenames <- factor(bac_otu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))
bac_zotu_elev <- bac_zotu
bac_zotu$Sitenames <- factor(bac_zotu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))

#make figure function based on elevation order instead
my.gg.elev <- function(df, yname,nudgenumber) {
  ggplot(data=df, aes(x=Sitenames, y=mean, col=Sitenames, label=Tukeylabels)) + geom_point(size=2) +
    labs(x=" ", y=yname, col="Site") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme_bw()+ #change to black and white theme
    theme(strip.text.x = element_text(size = 12, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
          axis.text.y=element_text(size=10),  #make y axis tick sizes bigger
          axis.title.y=element_text(size=14), #make y axis label larger
          panel.spacing.x = unit(0,"null"))+ #removing spacing between panels
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_color_manual(values=c("red", "green", "orange","blue","purple"))+#add in manual colors for points/lines
    scale_x_discrete(labels= c("275","470","1280","1710","2240"))+ #this change the x axis tick labels to just hte first letter instead of the entire word
    geom_text(nudge_y=nudgenumber) + #this is for the Tukey labels, makes them black and nudges them up on y axis so they aren't directly on top of point
    guides(col=FALSE) #remove legend for color
  #scale_y_continuous(position = "right") #put y axis on right side
}

#make figures 
p5 <- my.gg.elev(bac_otu_elev, "Bacterial OTU Richness",40) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p5

p6 <- my.gg.elev(bac_zotu_elev, "Bacterial ESV Richness",55) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p6

p7 <- my.gg.elev(fungi_otu_elev, "Fungal OTU Richness",10) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p7

p8 <- my.gg.elev(fungi_zotu_elev, "Fungal ESV Richness",15) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p8


#arrange just bacteria 
bacteria_otuvesv_elev <- ggarrange(p5,p6, ncol = 2, common.legend = TRUE, legend="bottom", labels=c("A","B"))
bacteria_otuvesv_elev

#arrange just fungi
fungi_otuvesv_elev <- ggarrange(p7,p8, ncol = 2, common.legend = TRUE, legend="bottom", labels=c("A","B")) #add titles A and B to plots
fungi_otuvesv_elev


pdf("Figures/zotuvotu_richness/Meanrichnessbysite_3timepoints_fungi_zotuvotu_elev.pdf", height=3, width=6)
fungi_otuvesv_elev
dev.off()

pdf("Figures/zotuvotu_richness/Meanrichnessbysite_3timepoints_bac_zotuvotu_elev.pdf", height=3, width=6)
bacteria_otuvesv_elev
dev.off()
