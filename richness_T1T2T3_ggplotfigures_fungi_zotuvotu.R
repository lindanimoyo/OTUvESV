#Nov 27, 2017
#updated Dec 18, 2017 to add in functions and make it more clean
#compare OTUs and ZOTUs

#Fungal richness 
#use same codes as fungi to make nice Figures and anovas
#Use alpha richness values from robert edgar usearch v10 commands, they are all rarefied to 10,000 sequences each
# Alpha diversity
#usearchv10 -alpha_div otutab_survey_e10000.txt -output alpha_otu_survey.txt
#usearchv10 -alpha_div zotutab_survey_e10000.txt -output alpha_zotu_survey.txt

#ITS2 library
#Reset R's Brain
rm(list=ls())

#install.packages("plyr")
library(plyr )
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")

###############################################################################################################
######load in OTU and zotu alpha reports - these are from usearch, before I've removed everything non-fungal though
#####################################################################################################################
#this is with the original bioinformatics, removing 77 bp
#check OTU v zotu comparisons
#T1
OTU_T1 <- read.csv("data/fungi_ITS2/alpha_otu_T1_fungi.csv", row.names=1)
zOTU_T1 <- read.csv("data/fungi_ITS2/alpha_zotuT1.csv", row.names=1)

head(zOTU_T1)
head(OTU_T1)

#row names matc
row.names(zOTU_T1) ==row.names(OTU_T1)

#remove everything that starts with Z to remove inoculum samples
zOTU_T1_2  <- zOTU_T1[-grep("^Z", row.names(zOTU_T1)), ]
OTU_T1_2  <- OTU_T1[-grep("^Z", row.names(OTU_T1)), ]

#row names matc
row.names(zOTU_T1_2) ==row.names(OTU_T1_2)


#remove everything that starts with Z to remove inoculum samples
zOTU_inoc  <- zOTU_T1[grep("^Z", row.names(zOTU_T1)), ]
OTU_inoc  <- OTU_T1[grep("^Z", row.names(OTU_T1)), ]

#row names matc
row.names(zOTU_inoc) ==row.names(OTU_inoc)

#T2
OTU_T2 <- read.csv("data/fungi_ITS2/alpha_otu_T2.csv", row.names=1)
zOTU_T2 <- read.csv("data/fungi_ITS2/alpha_zotuT2.csv", row.names=1)

head(zOTU_T2)
head(OTU_T2)

#row names matc
row.names(zOTU_T2) ==row.names(OTU_T2)

#T3
OTU_T3 <- read.csv("data/fungi_ITS2/alpha_otu_T3.csv", row.names=1)
zOTU_T3 <- read.csv("data/fungi_ITS2/alpha_zotuT3.csv", row.names=1)

head(zOTU_T3)
head(OTU_T3)

#row names matc
row.names(zOTU_T3) ==row.names(OTU_T3)

#survey
zOTU_survey <- read.csv("data/fungi_ITS2/alpha_zotu_survey1and2.csv", row.names=1)
OTU_survey <- read.csv("data/fungi_ITS2/alpha_otu_survey1and2.csv", row.names=1)
row.names(zOTU_survey) ==row.names(OTU_survey)

#remove everything that starts with Z to remove inoculum samples
zOTU_survey_2  <- zOTU_survey[-grep("^Z", row.names(zOTU_survey)), ]
OTU_survey_2  <- OTU_survey[-grep("^Z", row.names(OTU_survey)), ]

#row names matc
row.names(zOTU_survey_2) ==row.names(OTU_survey_2)


####################################################################################
#separate T1, T2, T3
####################################################################################
tablefunctionsurvey <- function(df){
  df$Type <- str_sub(row.names(df),1,2)
  df$Sample <- str_sub(row.names(df),4,4)
  df$DateSite <- str_sub(row.names(df),1,4)
  return(df)
}

tablefunctionransplant <- function(df){
  df$Type <- str_sub(row.names(df),1,2)
  df$Sample <- str_sub(row.names(df),1,1)
  return(df)
}

zOTU_survey_2 <- tablefunctionsurvey(zOTU_survey_2)
OTU_survey_2 <- tablefunctionsurvey(OTU_survey_2)

zOTU_inoc <- tablefunctionsurvey(zOTU_inoc)
OTU_inoc <- tablefunctionsurvey(OTU_inoc)

OTU_T1_2 <- tablefunctionransplant(OTU_T1_2)
zOTU_T1_2 <- tablefunctionransplant(zOTU_T1_2)

OTU_T2 <- tablefunctionransplant(OTU_T2)
zOTU_T2 <- tablefunctionransplant(zOTU_T2)

OTU_T3 <- tablefunctionransplant(OTU_T3)
zOTU_T3  <- tablefunctionransplant(zOTU_T3)

####################################################################################
#Get fungi mean sd and se with controls
####################################################################################
#get mean, SD, SE of T1 and T2 by site by inoculum
# Calculate the means, sd, n, and se.
averagefunction <- function(df, summarizeby){
  ddply(df, summarizeby, summarise,
        mean = mean(richness,na.rm=TRUE),
        sd = sd(richness, na.rm=TRUE),
        n = sum(!is.na( richness)),
        se = sd/sqrt(n)
  )
}


fungi_T1_zotu <- averagefunction(zOTU_T1_2,"Type")
fungi_T1_otu <- averagefunction(OTU_T1_2,"Type")

fungi_T2_zotu <- averagefunction(zOTU_T2,"Type")
fungi_T2_otu <- averagefunction(OTU_T2,"Type")

fungi_T3_zotu <- averagefunction(zOTU_T3,"Type")
fungi_T3_otu <- averagefunction(OTU_T3,"Type")

fungi_inoc_zotu <- averagefunction(zOTU_inoc,"Type")
fungi_inoc_otu <- averagefunction(OTU_inoc,"Type")

fungi_survey_zotu <- averagefunction(zOTU_survey_2,"DateSite")
fungi_survey_otu <- averagefunction(OTU_survey_2,"DateSite")

####################################################################################
#Make a dataframe for Figures with controls
####################################################################################
#####check that sample names are same for all 3
fungi_T1_otu$Type == fungi_T2_otu$Type

#combine rows
fungi_all_otu <- rbind(fungi_T1_otu,fungi_T2_otu,fungi_T3_otu,fungi_inoc_otu)
fungi_all_zotu <- rbind(fungi_T1_zotu,fungi_T2_zotu,fungi_T3_zotu,fungi_inoc_zotu)

#add site
fungi_all_otu$site <- str_sub(fungi_all_otu$Type,1,1)
fungi_all_otu$site[81:85] <-  c(1,2,3,5,4)
fungi_all_zotu$site <- str_sub(fungi_all_zotu$Type,1,1)
fungi_all_zotu$site[81:85] <-  c(1,2,3,5,4)

#add time points
fungi_all_otu$Timepoint <- c(rep("T1",nrow(fungi_T1_otu)),rep("T2",nrow(fungi_T2_otu)),rep("T3",nrow(fungi_T3_otu)),rep("Z",nrow(fungi_inoc_otu)))
fungi_all_zotu$Timepoint <- c(rep("T1",nrow(fungi_T1_zotu)),rep("T2",nrow(fungi_T2_zotu)),rep("T3",nrow(fungi_T3_zotu)),rep("Z",nrow(fungi_inoc_zotu)))

#make a list of Inoculum by substring from sample name
fungi_all_otu$Inoculum <- str_sub(fungi_all_otu$Type,2,2)
fungi_all_zotu$Inoculum <- str_sub(fungi_all_zotu$Type,2,2)

#remove controls
controls <- which(fungi_all_otu$Inoculum %in%c("C","L","N"))
fungi_all_otu <- fungi_all_otu[-controls, ]
fungi_all_zotu <- fungi_all_zotu[-controls, ]
#check datatable
fungi_all_otu
fungi_all_zotu

#make list of colors according to Jen's color scheme and add in some for controls
listofcolors1 <- c("red","green","blue","purple","orange")
#make list of colors that matches 
fungi_all_otu$colors <- c(rep(listofcolors1,15),listofcolors1)
fungi_all_zotu$colors <- c(rep(listofcolors1,15),listofcolors1)
#make ist of site names
sitenames1 <- c("Desert","Grassland","Pine-Oak","Subalpine","Scrubland")
sitenames2 <- c(rep("Desert",5),rep("Grassland",5),rep("Pine-Oak",5),rep("Scrubland",5),rep("Subalpine",5))
fungi_all_otu$sitenames <- c(rep(sitenames2,3), sitenames1)
fungi_all_zotu$sitenames <- c(rep(sitenames2,3), sitenames1)

#make factors for site names and inoculum in correct order so they show up correct on ggplot figure
factororderfunction <- function(df){
  df$site<- factor(df$site,levels=c(1,4,2,3,5))
  df$sitenames <- factor(df$sitenames,levels=c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))
  df$Inoculum <- factor(df$Inoculum,levels=c("D","W","G","P","S"))
  return(df)
}

fungi_all_otu <-factororderfunction(fungi_all_otu)
fungi_all_zotu <-factororderfunction(fungi_all_zotu)


##do same for survey
#change name of survey
names(fungi_survey_otu)[1] <- "Type"
names(fungi_survey_zotu)[1] <- "Type"
#survey otu
fungi_survey_otu$site <- str_sub(fungi_survey_otu$Type,4,4)
fungi_survey_otu$Timepoint <- str_sub(fungi_survey_otu$Type,1,2)
fungi_surveyinoc_otu <- rbind(fungi_all_otu[fungi_all_otu$Timepoint=="Z", 1:7], fungi_survey_otu)
#survey zotu
fungi_survey_zotu$site <- str_sub(fungi_survey_zotu$Type,4,4)
fungi_survey_zotu$Timepoint <- str_sub(fungi_survey_zotu$Type,1,2)
fungi_surveyinoc_zotu <- rbind(fungi_all_zotu[fungi_all_zotu$Timepoint=="Z", 1:7], fungi_survey_zotu)
#put factors in right order otu

factororderfunctionsurvey <- function(df){
  df$site<- factor(df$site,levels=c(6,1,4,2,3,5))
  df$Timepoint <- factor(df$Timepoint,levels=c("Z","T0","T1","T2","T3"))
  return(df)
}

fungi_surveyinoc_otu <- factororderfunctionsurvey(fungi_surveyinoc_otu )
fungi_surveyinoc_zotu <- factororderfunctionsurvey(fungi_surveyinoc_zotu )

inocs <- which(fungi_all_zotu$Timepoint=="Z")
fungi_all_zotu <-fungi_all_zotu[-inocs,]
fungi_all_otu <-fungi_all_otu[-inocs,]


####################################################################################
#export to csv
####################################################################################
write.csv(fungi_inoc_otu,"data/fungal_mean_inoculum_OTU.csv")

write.csv(fungi_surveyinoc_otu,"data/fungalrichness_surveyinoc_e10000_otu.csv")
write.csv(fungi_all_otu,"data/fungalrichness_transplant_e10000_otu.csv")

write.csv(fungi_surveyinoc_zotu,"data/fungalrichness_surveyinoc_e10000_zotu.csv")
write.csv(fungi_all_zotu,"data/fungalrichness_transplant_e10000_zotu.csv")

####################################################################################
#get max and min for otu v zotu
####################################################################################
#get mean, SD, SE of T1 and T2 by site by inoculum
max(fungi_all_zotu$mean) #183.75
min(fungi_all_zotu$mean) #73.5

max(fungi_all_otu$mean) #157.75
min(fungi_all_otu$mean) #55.5

mean(fungi_all_otu$mean) #95

####################################################################################
#Make ggplot Figures/richness/, for sites*inoculum, for all 3 time points face wrap with and without controls
####################################################################################
#make vectors of colors and names
colors1 <-c("red", "orange", "green","blue","purple")
colors2 <- c("red", "orange", "green","blue","purple","grey")
sitenameslabs <- c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine")
sitenameslabs2 <- c("Salton","Desert","Scrubland","Grassland","Pine-Oak","Subalpine")

####make figure function
ggplotrichnessfunction <- function(df, ylab, xvariable, colorvariable, groupvariable, xscalelabels,xscalecolors){
  ggplot(data=df, aes(x=xvariable, y=mean, col= colorvariable, group=groupvariable)) + geom_point(size=2) + geom_line()+ #group and geom_line add in the connector lines
    theme_bw() + #make black and white
    labs(x=" ", y=ylab, col="Inoculum") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
          axis.text.y=element_text(size=12),  #make y axis tick sizes bigger
          axis.title=element_text(size=14))+ #make y axis label larger  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_x_discrete(labels= xscalelabels)+
    scale_color_manual(labels=xscalelabels,#manual labels for legend
                       values=xscalecolors)   #add in manual colors for points/lines 
}

####change labels from T1 T2 T3 to 6, 12, 18 months
#make a labeller function a la https://stackoverflow.com/questions/3472980/ggplot-how-to-change-facet-labels
timepointnames <- c('T1'="6 months",'T2'="12 months",'T3'="18 months")

#run functions
p1_otu <- ggplotrichnessfunction(fungi_all_otu, "Fungal OTU richness", fungi_all_otu$site, fungi_all_otu$Inoculum,  fungi_all_otu$Inoculum,sitenameslabs,colors1) + facet_wrap(~Timepoint, ncol=3, labeller=as_labeller(timepointnames))
p1_zotu <- ggplotrichnessfunction(fungi_all_zotu, "Fungal ESV richness",fungi_all_zotu$site, fungi_all_zotu$Inoculum,  fungi_all_zotu$Inoculum,sitenameslabs,colors1) + facet_wrap(~Timepoint, ncol=3, labeller=as_labeller(timepointnames))

#plot and save OTU function
p1_otu


pdf("Figures/richness/fungalrichness_T1T2T2T3_otu.pdf", height=4, width=6)  
p1_otu 
dev.off()

#plot and save zOTU function
p1_zotu

pdf("Figures/richness/fungalrichness_T1T2T2T3_zotu.pdf", height=4, width=6)  
p1_zotu 
dev.off()

#plot and save otu v zotu
ggarrange(p1_otu,p1_zotu, ncol = 2, nrow = 1, common.legend = TRUE)
ggarrange(p1_otu,p1_zotu, ncol = 2, nrow = 1, common.legend = TRUE)%>%
  ggexport(filename = "Figures/richness/fungalrichness_T1T2T2T3_zotuvotu.pdf")


####################################################################################
#Now make Figures/richness/ with inoculum against time as x axis, and facet wrap the sites 
####################################################################################

#make ggplot of the mean % mass loss by site over time
timeplotfunction <- function(df, ylab){
  ggplot(data=df, aes(x=Timepoint, y=mean, col=Inoculum, group=Inoculum)) + geom_point(size=2) + geom_line()+ #group and geom_line add in the connector lines
    theme_bw() + #make black and white
    labs(x=" ", y=ylab, col="Inoculum") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
          axis.text.y=element_text(size=12),  #make y axis tick sizes bigger
          axis.title=element_text(size=14))+ #make y axis label larger
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_color_manual(labels=c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"),
                       values=c("red", "orange", "green","blue","purple"))+ #add in manual colors for points/lines+
  scale_x_discrete(labels=c ("T1" = "6","T2" = "12", "T3"="18"))  #change x axis scale to be 0, 6, 12, 18months
}

t_otu <- timeplotfunction(fungi_all_otu,"Fungal OTU richness")+ facet_wrap(~sitenames, ncol=3)
t_zotu <- timeplotfunction(fungi_all_zotu,"Fungal ESV richness")+ facet_wrap(~sitenames, ncol=3)

#plot and save OTU function
t_otu
t_otu %>%
  ggexport(filename = "Figures/richness/fungalrichness_bytime_otu.pdf")

#plot and save zOTU function
t_zotu
t_zotu %>%
  ggexport(filename = "Figures/richness/fungalrichness_bytime_zotu.pdf")

#plot and save otu v zotu
ggarrange(t_otu,t_zotu, ncol = 2, nrow = 1, common.legend = TRUE)
ggarrange(t_otu,t_zotu, ncol = 2, nrow = 1, common.legend = TRUE)%>%
  ggexport(filename = "Figures/richness/fungalrichness_bytime_zotuvotu.pdf")

####################################################################################
####make figure with inoc and survey
####################################################################################
ggplotrichnessfunction2 <- function(df, ylab,groupvariable,shapevariable, xscalelabels){
  ggplot(data=df, aes(x=site, y=mean, col= Timepoint, group=groupvariable, shape=shapevariable)) + geom_point(size=2) + geom_line()+ #group and geom_line add in the connector lines
    theme_bw() + #make black and white
    labs(x=" ", y=ylab, col="Timepoint") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
          axis.text.y=element_text(size=12),  #make y axis tick sizes bigger
          axis.title=element_text(size=14))+ #make y axis label larger  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_x_discrete(labels= xscalelabels)
}


p2_otu <- ggplotrichnessfunction2(fungi_surveyinoc_otu, "Fungal OTU richness",fungi_surveyinoc_otu$Timepoint,fungi_surveyinoc_otu$Timepoint, sitenameslabs2) + facet_wrap(~Timepoint, ncol=3)
p2_zotu <- ggplotrichnessfunction2(fungi_surveyinoc_zotu, "Fungal ESV richness",fungi_surveyinoc_zotu$Timepoint,fungi_surveyinoc_otu$Timepoint, sitenameslabs2) + facet_wrap(~Timepoint, ncol=3)

#plot and save OTU function
p2_otu 
p2_otu  %>%
  ggexport(filename = "Figures/richness/Fungalrichness_surveyandinoc_bysite_otu.pdf")

#plot and save zOTU function
p2_zotu
p2_zotu %>%
  ggexport(filename = "Figures/richness/Fungalrichness_surveyandinoc_bysite_zotu.pdf")

surveyrichnessplot <- ggarrange(p2_otu,p2_zotu, common.legend=TRUE)

surveyrichnessplot 

surveyrichnessplot  %>%
  ggexport(filename = "Figures/richness/Fungalrichness_survey_bytimepoint_otuvzotu.pdf")

#bind togetehr otu and zotu
fungi_surveyinoc_otu$amplicon <- rep("OTU", nrow(fungi_surveyinoc_otu))
fungi_surveyinoc_zotu$amplicon <- rep("ESV", nrow(fungi_surveyinoc_zotu))
fungi_surveyinoc_both <- rbind(fungi_surveyinoc_otu,fungi_surveyinoc_zotu)
fungi_surveyinoc_both$amplicon <- as.factor(fungi_surveyinoc_both$amplicon)

surveyboth <- ggplotrichnessfunction2(fungi_surveyinoc_both, "Fungal richness", fungi_surveyinoc_both$amplicon,fungi_surveyinoc_both$amplicon, sitenameslabs2) 
surveyboth+ facet_wrap(~Timepoint, ncol=3) 

#make figure with them both together
surveyboth <- ggplot(data=fungi_surveyinoc_both, aes(x=site, y=mean, col= Timepoint, group=amplicon, shape=amplicon)) + geom_point(size=4) + geom_line()+ #group and geom_line add in the connector lines
  theme_bw() + #make black and white
  labs(x=" ", y="Fungal richness", col="Timepoint") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
  theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=12),  #make y axis tick sizes bigger
        axis.title=element_text(size=14))+ #make y axis label larger  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
  scale_x_discrete(labels= sitenameslabs2)

surveyboth2 <- surveyboth + facet_wrap(~Timepoint, ncol=3) 

surveyboth2

surveyboth2  %>%
  ggexport(filename = "Figures/richness/Fungalrichness_survey_bytimepoint_otuvzotu_samefigure.pdf")


####################################################################################
#Now make Inoculum figures
####################################################################################

fungi_inoc_otu <- fungi_surveyinoc_otu[which(fungi_surveyinoc_otu$Timepoint=="Z"), ]
fungi_inoc_zotu <- fungi_surveyinoc_zotu[which(fungi_surveyinoc_zotu$Timepoint=="Z"), ]

#create a vector of Tukey labels based on above tukey tests
fungi_inoc_zotu$Tukeylabels <- c("a,b","b","a","c","c")
fungi_inoc_otu$Tukeylabels <- c("a,b","b","a","c","c")

####make figure with just inoc
ggplot_inocfunction <- function(df, ylab, nudgefactor){
   ggplot(data=df, aes(x=site, y=mean, col=site,label=Tukeylabels)) + geom_point(size=2) + geom_line()+ #group and geom_line add in the connector lines
    theme_bw() + #make black and white
    labs(x=" ", y=ylab, col="Inoculum") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
          axis.text.y=element_text(size=12),  #make y axis tick sizes bigger
          axis.title=element_text(size=14))+ #make y axis label larger
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_color_manual(labels=sitenameslabs, #manual labels for legend
                       values=colors1) + #add in manual colors for points/lines
    scale_x_discrete(labels= sitenameslabs )+ #change labels for x axis
    geom_text(nudge_y = nudgefactor)
}

inoc1 <- ggplot_inocfunction(fungi_inoc_otu, "Fungal OTU richness",40)
inoc2 <- ggplot_inocfunction(fungi_inoc_zotu, "Fungal ESV richness",40)

#plot and save OTU
inoc1
inoc1 %>%
  ggexport(filename = "Figures/richness/Fungalinoculum_richness_otu.pdf")

#plot and save zOTU
inoc2
inoc2 %>%
  ggexport(filename = "Figures/richness/Fungalinoculum_richness_zotu.pdf")

#plot and save OTU v zOTU

ggarrange(inoc1,inoc2, ncol = 2, common.legend = TRUE) 

ggarrange(inoc1,inoc2, ncol = 2, common.legend = TRUE) %>%
  ggexport(filename = "Figures/richness/Fungalrichness_inoc_zotuvotu.pdf")

####################################################################################
#Set up T1, T2, T3 for ANOVAs
####################################################################################
###########T1
par(mfrow=c(1,1))
####################################################################################
#Set up T1, T2, T3 for ANOVAs
####################################################################################
#remove controls from dataframe

T3controls <- which(zOTU_T3$Type %in% c("1C","2C","3C","4C","5C"))
zOTU_T3 <- zOTU_T3[-T3controls, ]
OTU_T3 <- OTU_T3[-T3controls, ]

anovasetupfunction <-function(df){
  df$Inoculum <- str_sub(row.names(df),2,2)
  df$Site <- str_sub(row.names(df),1,1)
  df$Inoculum <- factor(df$Inoculum,levels=c("D","W","G","P","S"))
  df$Site <- factor(df$Site,levels=c("1","4","2","3","5"))
  return(df)
}
###########T1
OTU_T1_2 <- anovasetupfunction(OTU_T1_2)
zOTU_T1_2 <- anovasetupfunction(zOTU_T1_2)
###########T2
OTU_T2 <- anovasetupfunction(OTU_T2)
zOTU_T2 <- anovasetupfunction(zOTU_T2)
###########T3
OTU_T3 <- anovasetupfunction(OTU_T3)
zOTU_T3 <- anovasetupfunction(zOTU_T3)

####################################################################################
#ANOVAs: Normal type 1 ANOVA
####################################################################################
#Do ANOVA model to test effect of inoculum and site
modelinoc<-aov(richness~Type,  data=OTU_inoc, na.action=na.omit)
summary(modelinoc)
summary.lm(modelinoc)


modelfungi_T1 <-aov(richness~Site*Inoculum,  data=OTU_T1_2, na.action=na.omit)
summary(modelfungi_T1)
summary.lm(modelfungi_T1)

modelfungi_T2<-aov(richness~Site*Inoculum,  data=OTU_T2, na.action=na.omit)
summary(modelfungi_T2)
summary.lm(modelfungi_T2)

modelfungi_T3<-aov(richness~Site*Inoculum,  data=OTU_T3, na.action=na.omit)
summary(modelfungi_T3)
summary.lm(modelfungi_T3)


#combine all time points for total anova
OTU_T1_2$Timepoint <- rep("T1",nrow(OTU_T1_2))
OTU_T2$Timepoint <- rep("T2",nrow(OTU_T2))
OTU_T3$Timepoint <- rep("T3",nrow(OTU_T3))

#bind them all together
OTU_T1T2T3 <- rbind(OTU_T1_2,OTU_T2,OTU_T3)

#make timepint a  factor
OTU_T1T2T3$Timepoint <- as.factor(OTU_T1T2T3$Timepoint)

#set up model
modelfungi_T1T2T3 <-aov(richness~Site*Inoculum*Timepoint,  data=OTU_T1T2T3, na.action=na.omit)
summary(modelfungi_T1T2T3 )
summary.lm(modelfungi_T1T2T3 )

##capture output
capture.output(summary.lm(modelfungi_T1T2T3 ),file="results/Richness_fungi_sitebytimepoint_otu.doc")
capture.output(summary(modelinoc),file="results/Richness_fungi_sitebyinoc_INOCULUM_otu.doc")
capture.output(summary(modelfungi_T1),file="results/Richness_fungi_sitebyinoc_T1_otu.doc")
capture.output(summary(modelfungi_T2),file="results/Richness_fungi_sitebyinoc_T2_otu.doc")
capture.output(summary(modelfungi_T3),file="results/Richness_fungi_sitebyinoc_T3_otu.doc")

##Zotu
#Do ANOVA model to test effect of inoculum and site
modelinoc_zotu<-aov(richness~Type,  data=zOTU_inoc, na.action=na.omit)
summary(modelinoc_zotu)
summary.lm(modelinoc_zotu)

modelfungi_T1_zotu <-aov(richness~Site*Inoculum,  data=zOTU_T1_2, na.action=na.omit)
summary(modelfungi_T1_zotu)
summary.lm(modelfungi_T1_zotu)

modelfungi_T2_zotu<-aov(richness~Site*Inoculum,  data=zOTU_T2, na.action=na.omit)
summary(modelfungi_T2_zotu)
summary.lm(modelfungi_T2_zotu)

modelfungi_T3_zotu<-aov(richness~Site*Inoculum,  data=zOTU_T3, na.action=na.omit)
summary(modelfungi_T3_zotu)
summary.lm(modelfungi_T3_zotu)


##capture output
capture.output(summary(modelinoc_zotu),file="results/Richness_fungi_sitebyinoc_INOCULUM_zotu.doc")
capture.output(summary(modelfungi_T1_zotu),file="results/Richness_fungi_sitebyinoc_T1_zotu.doc")
capture.output(summary(modelfungi_T2_zotu),file="results/Richness_fungi_sitebyinoc_T2_zotu.doc")
capture.output(summary(modelfungi_T3_zotu),file="results/Richness_fungi_sitebyinoc_T3_zotu.doc")

####################################################################################
#ANOVAs: Steve's way
####################################################################################
library(nlme)
m.1 <- gls(richness~Site*Inoculum,data=zOTU_T1_2,na.action="na.omit")
Anova(m.1,type=3)

m.2 <- gls(richness~Site*Inoculum,data=zOTU_T2,na.action="na.omit")
Anova(m.2,type=3)

m.3 <- gls(richness~Site*Inoculum,data=zOTU_T3,na.action="na.omit")
Anova(m.3,type=3)


####################################################################################
#Tukey HSD post hoc tests for Inoculum and Site OTU and ZOTU
####################################################################################
library(multcomp)
OTU_inoc$Site <-as.factor(OTU_inoc$Type)
OTU_inoc$Site <- factor(OTU_inoc$Site, levels=c("ZD","ZW","ZG","ZP","ZS"))

zOTU_inoc$Site <-as.factor(zOTU_inoc$Type)
zOTU_inoc$Site <- factor(zOTU_inoc$Site, levels=c("ZD","ZW","ZG","ZP","ZS"))

multcompfunction <- function(df, factorname){
   model <-aov(richness~factorname,  data=df)
   tuk <- glht(model, linfct = mcp(factorname = "Tukey")) 
   cld(tuk) 
 }

#by site: OTU
tuk.cld.inoc.otu.site <-multcompfunction(OTU_inoc, OTU_inoc$Site)
tuk.cld.T1.otu.site <-multcompfunction(OTU_T1_2, OTU_T1_2$Site)
tuk.cld.T2.otu.site <-multcompfunction(OTU_T2, OTU_T2$Site)
tuk.cld.T3.otu.site <-multcompfunction(OTU_T3, OTU_T3$Site)
#by site: zOTU
tuk.cld.inoc.zotu.site <-multcompfunction(zOTU_inoc, zOTU_inoc$Site)
tuk.cld.T1.zotu.site <-multcompfunction(zOTU_T1_2, zOTU_T1_2$Site)
tuk.cld.T2.zotu.site <-multcompfunction(zOTU_T2, zOTU_T2$Site)
tuk.cld.T3.zotu.site <-multcompfunction(zOTU_T3, zOTU_T3$Site)

#by Inoculum: OTU
tuk.cld.T1.otu.inoc <-multcompfunction(OTU_T1_2, OTU_T1_2$Inoculum)
tuk.cld.T2.otu.inoc <-multcompfunction(OTU_T2, OTU_T2$Inoculum)
tuk.cld.T3.otu.inoc <-multcompfunction(OTU_T3, OTU_T3$Inoculum)

#by Inoculum: zOTU
tuk.cld.T1.zotu.inoc <-multcompfunction(zOTU_T1_2, zOTU_T1_2$Inoculum)
tuk.cld.T2.zotu.inoc <-multcompfunction(zOTU_T2, zOTU_T2$Inoculum)
tuk.cld.T3.zotu.inoc <-multcompfunction(zOTU_T3, zOTU_T3$Inoculum)

#make figure for Tukey Site:OTU 
### make figure for OTU: use sufficiently large upper margin 
pdf("Figures/richness/Fungalrichness_T1T2T3_tukey_site_otu.pdf", height=7, width=8,pointsize=12)
old.par <- par(mai=c(1.1,0.8,1.25,0.1),mfrow=c(2,2),no.readonly = TRUE) #make enough space at top for tukey symbols and smaller spaces between Figures/richness/
plot(tuk.cld.T1.otu.site, ylab="Fungal OTU richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T1", line=5) #add in T1 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
plot(tuk.cld.T2.otu.site,  ylab="Fungal OTU richness",xaxt="n",xlab="",ylim=c(0,200)) 
mtext(side=3, "T2",line=5)#add in T2 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2)
plot(tuk.cld.T3.otu.site, ylab="Fungal OTU  richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T3", line=5) #add in T3 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
plot(tuk.cld.inoc.otu.site, ylab="Fungal OTU richness", xaxt="n",xlab="" ,ylim=c(0,550)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "Inoculum", line=5) #add in T3 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular

par(old.par)
dev.off()

#make figure for Tukey Site:zOTU 
### make figure for OTU: use sufficiently large upper margin 
pdf("Figures/richness/Fungalrichness_T1T2T3_tukey_site_zotu.pdf", height=7, width=8,pointsize=12)
old.par <- par(mai=c(1.1,0.8,1.25,0.1),mfrow=c(2,2),no.readonly = TRUE) #make enough space at top for tukey symbols and smaller spaces between Figures/richness/
plot(tuk.cld.T1.zotu.site, ylab="Fungal ESV richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T1", line=5) #add in T1 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
plot(tuk.cld.T2.zotu.site,  ylab="Fungal ESV richness",xaxt="n",xlab="",ylim=c(0,200)) 
mtext(side=3, "T2",line=5)#add in T2 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2)
plot(tuk.cld.T3.zotu.site, ylab="Fungal ESV richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T3", line=5) #add in T3 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
plot(tuk.cld.inoc.zotu.site, ylab="Fungal ESV richness", xaxt="n",xlab="" ,ylim=c(0,550)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "Inoculum", line=5) #add in T3 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular

par(old.par)
dev.off()


#make figure for Tukey Inoculum:OTU 
### make figure for OTU: use sufficiently large upper margin 
pdf("Figures/richness/Fungalrichness_T1T2T3_tukey_inoc_otu.pdf", height=7, width=8,pointsize=12)
old.par <- par(mai=c(1.1,0.8,1.25,0.1),mfrow=c(1,3),no.readonly = TRUE) #make enough space at top for tukey symbols and smaller spaces between Figures/richness/
plot(tuk.cld.T1.otu.inoc, ylab="Fungal OTU richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T1", line=5) #add in T1 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
plot(tuk.cld.T2.otu.inoc,  ylab="Fungal OTU richness",xaxt="n",xlab="",ylim=c(0,200)) 
mtext(side=3, "T2",line=5)#add in T2 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2)
plot(tuk.cld.T3.otu.inoc, ylab="Fungal OTU richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T3", line=5) #add in T3 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
par(old.par)
dev.off()

#make figure for Tukey Inoculum:zOTU
pdf("Figures/richness/Fungalrichness_T1T2T3_tukey_inoc_zotu.pdf", height=7, width=8,pointsize=12)
old.par <- par(mai=c(1.1,0.8,1.25,0.1),mfrow=c(1,3),no.readonly = TRUE) #make enough space at top for tukey symbols and smaller spaces between Figures/richness/
plot(tuk.cld.T1.zotu.inoc, ylab="Fungal ESV richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T1", line=5) #add in T1 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
plot(tuk.cld.T2.zotu.inoc,  ylab="Fungal ESV richness",xaxt="n",xlab="",ylim=c(0,200)) 
mtext(side=3, "T2",line=5)#add in T2 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2)
plot(tuk.cld.T3.zotu.inoc, ylab="Fungal ESV richness", xaxt="n",xlab="" ,ylim=c(0,200)) #add in same y limits for all, supress x axis labela nd tick marks
mtext(side=3, "T3", line=5) #add in T3 label
axis(side=1, at=c(1,2,3,4,5), labels=sitenameslabs, las=2) #add in customized x axis labels and make them perpendicular
par(old.par)
dev.off()



###################################################################################
#calculating effect sizes with eta squared: OTU
####################################################################################
#https://egret.psychol.cam.ac.uk/statistics/local_copies_of_sources_Cardinal_and_Aitken_ANOVA/glm_effectsize.htm
#https://artax.karlin.mff.cuni.cz/r-help/library/lsr/html/etaSquared.html

#install.packages("lsr")
library(lsr)

etasquaredT1T2T3 <- cbind(etaSquared(modelfungi_T1, type=2),etaSquared(modelfungi_T2, type=2),etaSquared(modelfungi_T3, type=2))

colnames(etasquaredT1T2T3) <- c("T1 eta.sq","T1 eta.sq.part", "T2 eta.sq","T2 eta.sq.part","T3 eta.sq","T3 eta.sq.part")

etasquaredT1T2T3
etasquaredT1T2T3trans <- t(etasquaredT1T2T3)
####################################################################################
#calculating effect sizes with omegasquared: OTU
####################################################################################
#source in functions
#https://gist.github.com/arnoud999/e677516ed45e9a11817e
source('~/Dropbox/StatsandProgramming/source/omegas.R', chdir = TRUE)

# Omega-squared using arnaud platinga code #https://gist.github.com/arnoud999/e677516ed45e9a11817e
Omegas(modelfungi_T2)
partialOmegas(modelfungi_T2)

#using code from here: https://stats.stackexchange.com/questions/2962/omega-squared-for-measure-of-effect-in-r
omega_sq(modelfungi_T2)

#all codes come out the exact same as Steve's except steve's has an error in it bc one of his come's out neg

#ok so now caluculate omegas for all 3
omegaT1 <- rbind(Omegas(modelfungi_T1),partialOmegas(modelfungi_T1))
row.names(omegaT1) <- c("omegasT1","partialomegasT1")

omegaT2 <- rbind(Omegas(modelfungi_T2),partialOmegas(modelfungi_T2))
row.names(omegaT2) <- c("omegasT2","partialomegasT2")

omegaT3 <- rbind(Omegas(modelfungi_T3),partialOmegas(modelfungi_T3))
row.names(omegaT3) <- c("omegasT3","partialomegasT3")

#combine all into one
omegasT1T2T3 <- rbind(omegaT1,omegaT2,omegaT3)
omegasT1T2T3

#combine with eta squared
omegasandetas <- rbind(omegasT1T2T3,etasquaredT1T2T3trans) 
omegasandetas

write.csv(omegasandetas, "results/fungi_richness_omegasandetasT1T2T3_otu.csv")


###################################################################################
#calculating effect sizes with eta squared: ZOTU
####################################################################################

etasquaredT1T2T3_zotu <- cbind(etaSquared(modelfungi_T1_zotu, type=2),etaSquared(modelfungi_T2_zotu, type=2),etaSquared(modelfungi_T3_zotu, type=2))

colnames(etasquaredT1T2T3_zotu) <- c("T1 eta.sq","T1 eta.sq.part", "T2 eta.sq","T2 eta.sq.part","T3 eta.sq","T3 eta.sq.part")

etasquaredT1T2T3_zotu
etasquaredT1T2T3_zotu_trans <- t(etasquaredT1T2T3_zotu)

#using code from here: https://stats.stackexchange.com/questions/2962/omega-squared-for-measure-of-effect-in-r

#ok so now caluculate omegas for all 3
omegaT1_zotu <- rbind(Omegas(modelfungi_T1_zotu),partialOmegas(modelfungi_T1_zotu))
row.names(omegaT1_zotu) <- c("omegasT1","partialomegasT1")

omegaT2_zotu <- rbind(Omegas(modelfungi_T2_zotu),partialOmegas(modelfungi_T2_zotu))
row.names(omegaT2_zotu) <- c("omegasT2","partialomegasT2")

omegaT3_zotu <- rbind(Omegas(modelfungi_T3_zotu),partialOmegas(modelfungi_T3_zotu))
row.names(omegaT3_zotu) <- c("omegasT3","partialomegasT3")

#combine all into one
omegasT1T2T3_zotu <- rbind(omegaT1_zotu,omegaT2_zotu,omegaT3_zotu)
omegasT1T2T3_zotu

#combine with eta squared
omegasandetas_zotu <- rbind(omegasT1T2T3_zotu,etasquaredT1T2T3_zotu_trans) 
omegasandetas_zotu

write.csv(omegasandetas_zotu, "results/fungi_richness_omegasandetasT1T2T3_zotu.csv")

class(omegasandetas_zotu[,1])

#check correlation for site
cor.test((omegasandetas_zotu[,1]),(omegasandetas[,1])) #95%
#check correlation for site
cor.test((omegasandetas_zotu[,2]),(omegasandetas[,2])) #92%
#check correlation for site
cor.test((omegasandetas_zotu[,3]),(omegasandetas[,3])) #98%
####################################################################################
#Make overall site mean figure with Tukey letters on it
####################################################################################
#run function to get average by Site 
fungi_T1_site_zotu <- averagefunction(zOTU_T1_2,"Site")
fungi_T1_site_otu <- averagefunction(OTU_T1_2,"Site")

fungi_T2_site_zotu <- averagefunction(zOTU_T2,"Site")
fungi_T2_site_otu <- averagefunction(OTU_T2,"Site")

fungi_T3_site_zotu <- averagefunction(zOTU_T3,"Site")
fungi_T3_site_otu <- averagefunction(OTU_T3,"Site")

#combine them all
fungi_T1T2T3_site_all_otu <- rbind(fungi_T1_site_otu,fungi_T2_site_otu,fungi_T3_site_otu)
fungi_T1T2T3_site_all_otu$Timepoint <- c(rep("T1",nrow(fungi_T1_site_otu)),rep("T2",nrow(fungi_T2_site_otu)),rep("T3",nrow(fungi_T3_site_otu)))
fungi_T1T2T3_site_all_otu$site <- factor(fungi_T1T2T3_site_all_otu$Site ,levels=c(1,4,2,3,5))
#create a vector of Tukey labels based on above tukey tests
fungi_T1T2T3_site_all_otu$Tukeylabels <- c("a","a","a","a","a","a","b","a,b","b","b","a","a,b","b","b","a,b")

#combine them all
fungi_T1T2T3_site_all_zotu <- rbind(fungi_T1_site_zotu,fungi_T2_site_zotu,fungi_T3_site_zotu)
fungi_T1T2T3_site_all_zotu$Timepoint <- c(rep("T1",nrow(fungi_T1_site_zotu)),rep("T2",nrow(fungi_T2_site_zotu)),rep("T3",nrow(fungi_T3_site_zotu)))
fungi_T1T2T3_site_all_zotu$site <- factor(fungi_T1T2T3_site_all_zotu$Site ,levels=c(1,4,2,3,5))
#create a vector of Tukey labels based on above tukey tests
fungi_T1T2T3_site_all_zotu$Tukeylabels <- c("a","a","a","a","a","a","a,b","a,b","b","a,b","a","a,b","a,b","b","a,b")


####ggplot figure for mean richness with function from earlier for otu
p5 <- ggplot_inocfunction(fungi_T1T2T3_site_all_otu,"Fungal OTU Richness",10) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))

####ggplot figure for mean richness with function from earlier for zotu
p6 <- ggplot_inocfunction(fungi_T1T2T3_site_all_zotu,"Fungal ESV Richness",15) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))

range(fungi_T1T2T3_site_all_otu$mean)
mean(fungi_T1T2T3_site_all_otu$mean)
#look at plots and save output : OTU
p5
pdf("Figures/richness/MeanfungirichnessbySitebytimepoint_facetwrap_otu.pdf", height=4, width=6)
p5 
dev.off()

#look at plots and save output : xOTU
pdf("Figures/richness/MeanfungirichnessbySitebytimepoint_facetwrap_zotu.pdf", height=4, width=6)
p6 
dev.off()


#put them together:
ggarrange(p5,p6, ncol = 2, nrow = 1, common.legend = TRUE)
ggarrange(p5,p6, ncol = 2, nrow = 1, common.legend = TRUE) %>%
  ggexport(filename = "Figures/richness/MeanfungirichnessbySitebytimepoint_facetwrap_zotuvotu.pdf")
