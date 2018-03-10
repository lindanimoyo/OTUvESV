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
library(plyr )

#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")

################################################################################################
######Fungi: avg.dist bray-curtis dissimilarity matrices - median and square root transformed for T1 T2 T3 trasplant
################################################################################################
#bray curtis
fungi_T3_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T3_otu.csv", row.names=1 )
fungi_T3_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T3_zotu.csv", row.names=1 )

#check that row and column names all match
row.names(fungi_T3_bray_otu) == row.names(fungi_T3_bray_zotu)
row.names(fungi_T3_bray_otu) ==colnames(fungi_T3_bray_otu)
colnames(fungi_T3_bray_otu) ==colnames(fungi_T3_bray_zotu)

################################################################################################
######Bacteria: avg.dist bray-curtis dissimilarity matrices - median and square root transformed for T1 T2 T3 trasplant
################################################################################################
#bray-curtis
bac_T3_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T3_otu.csv", row.names=1 )
bac_T3_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T3_zotu.csv", row.names=1 )

#check that row and column names all match
row.names(bac_T3_bray_otu) == row.names(bac_T3_bray_zotu)
row.names(bac_T3_bray_otu) ==colnames(bac_T3_bray_otu)
colnames(bac_T3_bray_otu) ==colnames(bac_T3_bray_zotu)

######################################################################################################################
###### Make example NMDS plot and stress plots - Bac and Fungi, OTU and ZOTU
######################################################################################################################
#make the nmds 
plotnmds_fungi_T3_otu <- metaMDS(fungi_T3_bray_otu, k=2, trymax=100)
plotnmds_fungi_T3_otu #stress = 0.15
stressplot(plotnmds_fungi_T3_otu)

plotnmds_fungi_T3_zotu <- metaMDS(fungi_T3_bray_zotu, k=2, trymax=100)
plotnmds_fungi_T3_zotu #stress= 0.14
stressplot(plotnmds_fungi_T3_zotu)

plotnmds_bac_T3_otu <- metaMDS(bac_T3_bray_otu, k=3, trymax=100)
plotnmds_bac_T3_otu  #stress = 0.067
stressplot(plotnmds_bac_T3_otu)

plotnmds_bac_T3_zotu <- metaMDS(bac_T3_bray_zotu, k=3, trymax=100)
plotnmds_bac_T3_zotu #stress = 0.078
stressplot(plotnmds_bac_T3_zotu)

######################################################################################################################
###### read in predictor files and match thtem
######################################################################################################################

#read in predictor files
map_T3_bac <- read.csv("data/16S_usearch10/16S_T3_map.csv", row.names=1)

##Discard all samples in map that don't appear otu table
map_T3_bac<- map_T3_bac [map_T3_bac$samplename2 %in% rownames(bac_T3_bray_otu), ]
dim(map_T3_bac)
dim(bac_T3_bray_otu)

##Match row names in tables
map_T3_bac<- map_T3_bac[match(rownames(bac_T3_bray_otu),map_T3_bac$samplename2),]
map_T3_bac$samplename2 == rownames(bac_T3_bray_otu)
rownames(bac_T3_bray_otu) ==colnames(bac_T3_bray_otu)

#read in predictor files
map_T3_fungi <- read.csv("data/fungi_ITS2/map_T3_transplant.csv", row.names=1)

##Discard all samples in map that don't appear otu table
map_T3_fungi<- map_T3_fungi [map_T3_fungi$samplename2 %in% rownames(fungi_T3_bray_otu), ]
dim(map_T3_fungi)
dim(fungi_T3_bray_otu)

##Match row names in tables
map_T3_fungi<- map_T3_fungi[match(rownames(fungi_T3_bray_otu),map_T3_fungi$samplename2  ),]
map_T3_fungi$samplename2 == rownames(fungi_T3_bray_otu)
rownames(fungi_T3_bray_otu) ==colnames(fungi_T3_bray_otu)

######################################################################################################################
###### colored by site
######################################################################################################################
names(map_T3_fungi)
map_T3_fungi$colorbysite <- as.character(map_T3_fungi$colorbysite )
map_T3_fungi$colorbyinoc <- as.character(map_T3_fungi$colorbyinoc )

names(map_T3_bac)
map_T3_bac$colorbysite <- as.character(map_T3_bac$colorbysite )
map_T3_bac$colorbyinoc <- as.character(map_T3_bac$colorbyinoc )

plot(scores(plotnmds_fungi_T3_otu), col=map_T3_fungi$colorbyinoc, pch=map_T3_fungi$shapesbysite, cex=1)
plot(scores(plotnmds_fungi_T3_zotu), col=map_T3_fungi$colorbyinoc, pch=map_T3_fungi$shapesbysite, cex=1)

plot(scores(plotnmds_bac_T3_otu), col=map_T3_bac$colorbysite, pch=map_T3_bac$ShapesbyInoculum, cex=1)
plot(scores(plotnmds_bac_T3_zotu), col=map_T3_bac$colorbysite, pch=map_T3_bac$ShapesbyInoculum, cex=1)

######################################################################################################################
###### Figure for Publication A) colored by site B) colored by Inoc with bray curtis sqrt OTU
######################################################################################################################
plot.new()

#change it to make it simpler, no shapes
#make panels C and D
pdf("Figures/otuvzotu_dissimilarity/NMDS_ITS2_T3_OTUvESV.pdf",height=6, width=10)
par(mfrow=c(1,2), mai=c(1,0.8,1,0.2))
#OTU: plot the scores, color by inoculum, pch 16 ffor all, same y and x limits
plot(scores(plotnmds_fungi_T3_otu), col=map_T3_fungi$colorbyinoc, pch=16, cex=1, ylim=c(-0.4,0.6), xlim=c(-0.5,0.6))
#add in stress
mtext("Stress = 0.15", line=-2, adj=0.9, cex=1, side=3)
#add in OTU
mtext("Fungal OTU", line=1, adj=0.5, cex=2, side=3)
#add in panel letter and make it bold
mtext(expression(bold("C")), line=1, adj=0.05, cex=2, side=3)
#ESV: plot the scores, color by inoculum, pch 16 ffor all, same y and x limits
plot(scores(plotnmds_fungi_T3_zotu), col=map_T3_fungi$colorbyinoc, pch=16, cex=1,ylim=c(-0.4,0.6), xlim=c(-0.5,0.6))
#add in stress
mtext("Stress = 0.14", line=-2, adj=0.9, cex=1, side=3)
#add in OTU
mtext("Fungal ESV", line=1, adj=0.5, cex=2, side=3)
#add in panel letter and make it bold
mtext(expression(bold("D")), line=1, adj=0.05, cex=2, side=3)
dev.off()

pdf("Figures/otuvzotu_dissimilarity/NMDS_16S_T3_OTUvESV.pdf",height=6, width=10)
par(mfrow=c(1,2), mai=c(1,0.8,1,0.2))
#OTU: plot the scores, color by site, pch 16 ffor all, same y and x limits
plot(scores(plotnmds_bac_T3_otu), col=map_T3_bac$colorbysite, pch=16, cex=1, ylim=c(-0.4,0.4), xlim=c(-0.55,0.5))
#add in stress
mtext("Stress = 0.067", line=-2, adj=0.9, cex=1, side=3)
#add in OTU
mtext("Bacterial OTU", line=1, adj=0.5, cex=2, side=3)
#add in panel letter and make it bold
mtext(expression(bold("C")), line=1, adj=0.05, cex=2, side=3)
#ESV: plot the scores for the ZOTU, color by site, pch 16 for all, same x and y limits
plot(scores(plotnmds_bac_T3_zotu), col=map_T3_bac$colorbysite, pch=16, cex=1,ylim=c(-0.4,0.4),xlim=c(-0.55,0.5))
#add in stress
mtext("Stress = 0.078", line=-2, adj=0.9, cex=1, side=3)
#add in zotu
mtext("Bacterial ESV", line=1, adj=0.5, cex=2, side=3)
#add in panel letter and make it bold
mtext(expression(bold("D")), line=1, adj=0.05, cex=2, side=3)
dev.off()

par(mfrow=c(1,1),par(mai=c(1.02,0.82,0.82,0.42)))


######################################################################################################################
###### Make same figure with ggplot so I can arrange it with the richness figure 
######################################################################################################################
#make data frame for fungi 
fungi_T3_otu <- as.data.frame(scores(plotnmds_fungi_T3_otu))
fungi_T3_zotu <- as.data.frame(scores(plotnmds_fungi_T3_zotu))

#bind them togeter and get string for site and inoculum and put in correct order
fungi_T3_otuvzotu <- rbind(fungi_T3_otu,fungi_T3_zotu)

#make a column for Amplicon
fungi_T3_otuvzotu$Amplicon <- c(rep("OTU",nrow(fungi_T3_otu)), rep("ESV",nrow(fungi_T3_zotu)))
#order it to be OTU first
fungi_T3_otuvzotu$Amplicon <- factor(fungi_T3_otuvzotu$Amplicon,levels=c("OTU","ESV"))
#add timepoint
fungi_T3_otuvzotu$Timepoint <- rep("T3",nrow(fungi_T3_otuvzotu))
#add kingdom
fungi_T3_otuvzotu$Kingdom <- rep("Fungi",nrow(fungi_T3_otuvzotu))

#make function to get data frame with site and inoculum in correct order
dataframefunction <- function(df) {
  df$Site <- str_sub(row.names(df ), 2, 2)
  df$Inoculum <- str_sub(row.names(df), 3, 3)
  df$Site <- factor(df$Site,levels=c(1,4,2,3,5))
  df$Inoculum <- factor(df$Inoculum,levels=c("D","W","G","P","S"))
  return(df)
}

fungi_T3_otuvzotu <- dataframefunction(fungi_T3_otuvzotu)
fungi_T3_otuvzotu 


##############################
#make data frame for bacteria
##############################
bac_T3_otu <- as.data.frame(scores(plotnmds_bac_T3_otu))
bac_T3_zotu <- as.data.frame(scores(plotnmds_bac_T3_zotu))

#bind them togeter and get string for site and inoculum and put in correct order
bac_T3_otuvzotu <- rbind(bac_T3_otu,bac_T3_zotu)

#make a column for Amplicon
bac_T3_otuvzotu$Amplicon <- c(rep("OTU",nrow(bac_T3_otu)), rep("ESV",nrow(bac_T3_zotu)))
#order it to be OTU first
bac_T3_otuvzotu$Amplicon <- factor(bac_T3_otuvzotu$Amplicon,levels=c("OTU","ESV"))
#add timepoint
bac_T3_otuvzotu$Timepoint <- rep("T3",nrow(bac_T3_otuvzotu))
#add kingdom
bac_T3_otuvzotu$Kingdom <- rep("Bacteria",nrow(bac_T3_otuvzotu))

bac_T3_otuvzotu <- dataframefunction(bac_T3_otuvzotu)
bac_T3_otuvzotu 

#save files
write.csv(bac_T3_otuvzotu, "Figures/otuvzotu_dissimilarity/bacteria_NMDS_scores_otuzotu.csv")
write.csv(fungi_T3_otuvzotu , "Figures/otuvzotu_dissimilarity/fungi_NMDS_scores_otuzotu.csv")
######################################################################################################################
###### Make NMDS figure function in ggplot - colored by site and shapes by inoculum
######################################################################################################################
#make vectors for legend
sitenames <- c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine")
sitecolors  <- c("red","orange","green","blue","purple")
shapes <- c(17,3,16,15,18)

#includes an x axis
NMDS_sitecolor_function <- function(df, figurename) {
  #make ggplot object
  ggplot(df, aes(x=NMDS1, y=NMDS2, col=Inoculum, group=Inoculum)) +
    geom_point(size=2) +
    theme_bw() +
    labs(title = figurename)+
    scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
    #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
    theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text=element_text(size=10,angle=70, hjust=1),  #change size of  x and y axis tick labels
          axis.title=element_text(size=12),#change size of x and y axis labels
          legend.spacing = unit(0,"cm"), #change spacing between legends
          legend.text=element_text(size=10), #change size of legend text
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #remove major and minor grid panels
          legend.position="none") #remove legend
  
}

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
panelnamesfungi <- c('OTU'="C.                   Fungal OTU",'ESV'="D.                    Fungal ESV")

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
panelnamesbac <- c('OTU'="C.                   Bacterial OTU",'ESV'="D.                    Bacterial ESV")


#make NMDS for colors site and shapes inoculum
#bacteria all time points faceted
bacterialplot + facet_wrap(~Amplicon,labeller=as_labeller(panelnamesbac)) +stat_ellipse() #include ellipse around the group


