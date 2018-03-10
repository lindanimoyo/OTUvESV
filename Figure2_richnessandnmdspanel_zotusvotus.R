#Feb 19, 2018
#make beta diversity panels for fungi and bacteria separate OTU v ESV
#try it in base R and with ggplot

#Reset R's Brain
rm(list=ls())

#load libraries
library(plyr )
library(vegan)
library(ggplot2)
library(ggpubr)
library(stringr)
library(plyr )

#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")


###############################################################################################################
######First Panel: Richness Figure
###############################################################################################################
#read in fungi otu and zotu data frames with mean richness per site
fungi_otu <- read.csv("data/fungi_T1T2T3_meanrichnessbysite_otu.csv", row.names=1)
fungi_zotu <-  read.csv("data/fungi_T1T2T3_meanrichnessbysite_zotu.csv",  row.names=1)

#read in bacteria otu and zotu data frames with mean richness per site
bac_otu <- read.csv("data/bac_T1T2T3_site_all_otu.csv", row.names=1)
bac_zotu <- read.csv("data/bac_T1T2T3_meanrichnessbysite_zotu.csv",row.names=1)

#make factors correct order, based on elevation order
fungi_otu$Sitenames <- factor(fungi_otu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))
fungi_zotu$Sitenames <- factor(fungi_zotu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))
bac_otu$Sitenames <- factor(bac_otu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))
bac_zotu$Sitenames <- factor(bac_zotu$Sitenames ,levels=c("Desert","Grassland","Scrubland","Pine-Oak","Subalpine"))

###############################################################################################################
######Make figure function
##############################################################################################################

#includes an x axis
#make figure function based on elevation order instead
my.gg.elev <- function(df, yname,nudgenumber) {
  ggplot(data=df, aes(x=Sitenames, y=mean, col=Sitenames, label=Tukeylabels)) + geom_point(size=2) +
    labs(x=" ", y=yname, col="Site") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
    theme_bw()+ #change to black and white theme
    theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
          #axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
          axis.title.x = element_blank(),  #remove x axis title
          axis.text.x = element_blank(), #remove x axis text
          axis.ticks.x = element_blank(), #remove x axis ticks
          axis.text.y=element_text(size=10),  #make y axis tick sizes bigger
          axis.title.y=element_text(size=14), #make y axis label larger
          panel.spacing.x = unit(0,"null"))+ #removing spacing between panels
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
    scale_color_manual(values=c("red", "green", "orange","blue","purple"))+#add in manual colors for points/lines
    #scale_x_discrete(labels= c("275","470","1280","1710","2240"))+ #this change the x axis tick labels to elevations
    geom_text(nudge_y=nudgenumber) + #this is for the Tukey labels, makes them black and nudges them up on y axis so they aren't directly on top of point
    guides(col=FALSE) #remove legend for color
  #scale_y_continuous(position = "right") #put y axis on right side
}


####change labels from T1 T2 T3 to 6, 12, 18 months
#make a labeller function a la https://stackoverflow.com/questions/3472980/ggplot-how-to-change-facet-labels
timepointnames <- c('T1'= "6 mos",'T2'="12 mos",'T3'="18 mos")


#make figures 
p5 <- my.gg.elev(bac_otu, "Bacterial OTU Richness",40) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p5

p6 <- my.gg.elev(bac_zotu, "Bacterial ESV Richness",55) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p6

p7 <- my.gg.elev(fungi_otu, "Fungal OTU Richness",10) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p7

p8 <- my.gg.elev(fungi_zotu, "Fungal ESV Richness",15) + facet_wrap(~Timepoint, ncol=3,labeller=as_labeller(timepointnames))
p8


#arrange just bacteria 
bacteria_otuvesv_elev <- ggarrange(p5,p6, ncol = 2, common.legend = TRUE, legend="bottom")
bacteria_otuvesv_elev

#arrange just fungi
fungi_otuvesv_elev <- ggarrange(p7,p8, ncol = 2, common.legend = TRUE, legend="bottom") #add titles A and B to plots
fungi_otuvesv_elev


###############################################################################################################################################
######Read in NMDS scores files for ggplot for bacteria and fungi, rarefied and square root transformed bray-curtis dissimilarity matrices
############################################################################################################
#these are the data frames of teh NMDS scores from metaDS on the square root transformed rarefied bray-curtis dissimilarity matrices
#from this script: NMDS_zotusvotus_16SandITS2_T3.R
bac_T3_otuvzotu <- read.csv("Figures/otuvzotu_dissimilarity/bacteria_NMDS_scores_otuzotu.csv",row.names=1 )
fungi_T3_otuvzotu <- read.csv("Figures/otuvzotu_dissimilarity/fungi_NMDS_scores_otuzotu.csv",row.names=1 )

#stress for fungi OTU = 0.15
#stress for fungi ESV = 0.14
#stress for bacteria OTU = 0.067
#stress for bacteria ESV=0.078

#make function to get data frame with site and inoculum in correct order
#order correctly
bac_T3_otuvzotu$Amplicon <- factor(bac_T3_otuvzotu$Amplicon,levels=c("OTU","ESV"))
bac_T3_otuvzotu$Site <- factor(bac_T3_otuvzotu$Site,levels=c(1,4,2,3,5))
bac_T3_otuvzotu$Inoculum  <- factor(bac_T3_otuvzotu$Inoculum,levels=c("D","W","G","P","S"))
fungi_T3_otuvzotu$Amplicon <- factor(fungi_T3_otuvzotu$Amplicon,levels=c("OTU","ESV"))
fungi_T3_otuvzotu$Site  <- factor(fungi_T3_otuvzotu$Site,levels=c(1,4,2,3,5))
fungi_T3_otuvzotu$Inoculum <- factor(fungi_T3_otuvzotu$Inoculum,levels=c("D","W","G","P","S"))
######################################################################################################################
######GGplot Facet wrapping: Make NMDS figure function in ggplot - colored by site and shapes by inoculum
######################################################################################################################
#make vectors for legend
sitenames <- c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine")
sitecolors  <- c("red","orange","green","blue","purple")
shapes <- c(17,3,16,15,18)

fungalplot <- ggplot(fungi_T3_otuvzotu, aes(x=NMDS1, y=NMDS2, col=Inoculum, group=Inoculum)) + geom_point(size=2) +
  theme_bw() +
  labs(title = "Fungal NMDS")+
  scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=12),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #remove major and minor grid panels
        legend.position="none") #remove legend

fungalplot + facet_wrap(~Amplicon) +stat_ellipse() #include ellipse around the group

#make panel names function
panelnamesfungi <- c('OTU'="Fungal OTU",'ESV'="Fungal ESV")

fungalplot + facet_wrap(~Amplicon,labeller=as_labeller(panelnamesfungi )) +stat_ellipse() #include ellipse around the group


#make bacteria plot
bacterialplot <- ggplot(bac_T3_otuvzotu, aes(x=NMDS1, y=NMDS2, col=Site, group=Site)) + geom_point(size=2) +
  theme_bw() +
  #labs(title = c("Bacterial NMDS")+
  scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=12),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #remove major and minor grid panels
        legend.position="none") #remove legend

bacterialplot + facet_wrap(~Amplicon) +stat_ellipse() #include ellipse around the group


#make panel names function
panelnamesbac <- c('OTU'="Bacterial OTU",'ESV'="Bacterial ESV")


#make NMDS for colors site and shapes inoculum
#bacteria all time points faceted
bacterialplot + facet_wrap(~Amplicon,labeller=as_labeller(panelnamesbac)) +stat_ellipse() #include ellipse around the group


#save these as times
fungiNMDS <- fungalplot + facet_wrap(~Amplicon,labeller=as_labeller(panelnamesfungi )) +stat_ellipse() #include ellipse around the group
fungiNMDS

bacNMDS <- bacterialplot + facet_wrap(~Amplicon,labeller=as_labeller(panelnamesbac)) +stat_ellipse() #include ellipse around the group
bacNMDS
######################################################################################################################
###### Instead of facet wrapping make them individual: Bacteria
######################################################################################################################
#separate out OTU v ESV
bac_T3_otuvzotu_1 <- bac_T3_otuvzotu[which(bac_T3_otuvzotu$Amplicon=="OTU"), ]
bac_T3_otuvzotu_2 <- bac_T3_otuvzotu[which(bac_T3_otuvzotu$Amplicon=="ESV"), ]


bacterialplot1 <- ggplot(bac_T3_otuvzotu_1, aes(x=NMDS1, y=NMDS2, col=Site, group=Site)) + geom_point(size=2) +
  theme_bw() + #make black and white
  ylim(-0.5,0.5) + #set ylimits
  stat_ellipse()+  #add in ellipses
  geom_text(label="Stress = 0.067",y=0.5, x=0.3, col="black", size=3)+
 # labs(title = c("Bacterial OTU Bray-Curtis Dissimilarity at 18 mos; Stress = 0.067"))+
  labs(title = c("Bacterial OTU Bray-Curtis Dissimilarity"))+
  scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  theme(axis.text=element_text(size=12),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis label
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #remove major and minor grid panels
        legend.position="none") #remove legend


bacterialplot1 

bacterialplot2 <- ggplot(bac_T3_otuvzotu_2, aes(x=NMDS1, y=NMDS2, col=Site, group=Site)) + geom_point(size=2) +
  theme_bw() + #make black and white
  ylim(-0.5,0.5) + #set ylimits
  stat_ellipse()+  #add in ellipses
  geom_text(label="Stress = 0.078",y=0.5, x=0.3, col="black", size=3)+
  #labs(title = c("Bacterial ESV Bray-Curtis Dissimilarity at 18 mos; Stress = 0.078"))+
  labs(title = c("Bacterial ESV Bray-Curtis Dissimilarity"))+
  scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  theme(axis.text=element_text(size=12),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis label
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #remove major and minor grid panels
        legend.position="none") #remove legend

bacterialplot2

######################################################################################################################
###### Instead of facet wrapping make them individual: Fungi
#####################################################################################################################

fungi_T3_otuvzotu_1 <- fungi_T3_otuvzotu[which(fungi_T3_otuvzotu$Amplicon=="OTU"), ]
fungi_T3_otuvzotu_2 <- fungi_T3_otuvzotu[which(fungi_T3_otuvzotu$Amplicon=="ESV"), ]


fungalplot1 <- ggplot(fungi_T3_otuvzotu_1, aes(x=NMDS1, y=NMDS2, col=Inoculum, group=Inoculum)) + geom_point(size=2) +
  theme_bw() + #make black and white
  stat_ellipse()+  #add in ellipses
  ylim(-0.55,0.5)+#set ylimits
  geom_text(label="Stress = 0.15",y=0.5, x=0.6, col="black", size=3)+
  #labs(title = c("Fungal OTU Bray-Curtis Dissimilarity at 18 mos; Stress = 0.15"))+
  labs(title = c("Fungal OTU Bray-Curtis Dissimilarity"))+
  scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  theme(axis.text=element_text(size=12),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis label
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #remove major and minor grid panels
        legend.position="none") #remove legend


fungalplot1 

fungalplot2 <- ggplot(fungi_T3_otuvzotu_2, aes(x=NMDS1, y=NMDS2, col=Inoculum, group=Inoculum)) + geom_point(size=2) +
  theme_bw() + #make black and white
  stat_ellipse()+  #add in ellipses
  geom_text(label="Stress = 0.14",y=0.5, x=0.6, col="black", size=3)+
  ylim(-0.55,0.5)+ #set ylimits
  #labs(title = c("Fungal ESV Bray-Curtis Dissimilarity at 18 mos; Stress = 0.14"))+
  labs(title = c("Fungal ESV Bray-Curtis Dissimilarity"))+
  scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  theme(axis.text=element_text(size=12),  #change size of  x and y axis tick labels
        axis.title=element_text(size=14),#change size of x and y axis label
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #remove major and minor grid panels
        legend.position="none") #remove legend


fungalplot2
######################################################################################################################
######  Arrange ggplot items for a single figure
######################################################################################################################

bacteriapanel <- ggarrange(p5,p6,bacterialplot1,bacterialplot2,ncol=2,nrow=2, labels=c("A","B","C","D"))
bacteriapanel

pdf("Figures/zotuvotu_richness/Figure2_bacteria_richnessandNMDS_otuvzotu.pdf", height=6, width=8)
bacteriapanel
dev.off()

fungalpanel <- ggarrange(p7,p8,fungalplot1,fungalplot2,ncol=2,nrow=2, labels=c("A","B","C","D"))
fungalpanel 

pdf("Figures/zotuvotu_richness/Figure2_fungi_richnessandNMDS_otuvzotu.pdf", height=6, width=8)
fungalpanel 
dev.off()