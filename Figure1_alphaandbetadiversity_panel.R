
#Jan 30, 2018

# Make Figure 1 for otu v esv commentary
#combine alpha and beta diversity into one figure

#make figure 1 for OTU v ZOTU commentary
#panel A bacteria OTU v ZOTU one example correlation
#Panel B fungi OTU v ZOTU
#panel C boxplots/histograms of all corelations
#Reset R's Brain
rm(list=ls())

#install.packages("plyr")
library(plyr )
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(vegan)
#set working directory
setwd("~/Dropbox/StatsandProgramming/16SElevationGradient/")

bacrichness <- read.csv( "data/16S_otuvzotu.csv", row.names=1)
fungalrichness <- read.csv("data/ITS2_OTUvZOTUall.csv", row.names=1)

###############################################################################################################
######Make just bacteria figure
#####################################################################################################################

#change levels to 6, 12, 18 months
levels(bacrichness$Timepoint) <- c("6 months","12 months","18 months")

#make figure
bacscatter <- ggscatter(bacrichness, x = "richness", y = "zOTU_richness",
                        add = "reg.line",                         # Add regression line
                        conf.int = TRUE,                          # Add confidence interval
                        color = "Timepoint",            # Color by groups "cyl"
                        shape = "Timepoint",                             # Change point shape by groups "cyl"
                        size=2
)+
  #stat_cor(aes(color = Timepoint), label.x = 200, label.y.npc=0.1, method="pearson")     + # Add correlation coefficient
  labs(x="Bacterial OTU alpha-diversity",y="Bacterial ESV alpha-diversity", title="A") + # change labels
  theme(axis.text.x=element_text(size=14, angle=70, hjust=1), axis.text.y=element_text(size=14), #make y axis tick sizes bigger
        axis.title=element_text(size=14),
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger ) #make y axis label larger

bacscatter

bacscatter2 <-  bacscatter +   theme(legend.position="none")#remove legend
###############################################################################################################
######Make just fungifigure
#####################################################################################################################

#figure out what colors ggplot2 used to build this figure
ggplot_build(richness_fungi)$data[[1]][[1]]

#useful webpage: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/#comments-list

#change the levels from T1 T2 T3 to 6,12,18 months
levels(fungalrichness$Timepoint) <- c("6 months","12 months","18 months")

fungalscatter <- ggscatter(fungalrichness, x = "richness", y = "zOTU_richness",
                           add = "reg.line",                         # Add regression line
                           conf.int = TRUE,                          # Add confidence interval
                           color = "Timepoint",            # Color by groups "cyl"
                           shape = "Timepoint",                            # Change point shape by groups "cyl"
                           size=2
)+
  #stat_cor(aes(color = Timepoint),label.x = 100,label.y.npc=0.1, method="pearson")   + # Add correlation coefficient
  labs(x="Fungal OTU alpha-diversity",y="Fungal ESV alpha-diversity", title="B") + # change labels
  theme(axis.title=element_text(size=14), axis.text.x=element_text(size=14, angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger  

fungalscatter 

fungalscatter2 <- fungalscatter +   theme(legend.position="none")#remove legend

###############################################################################################################
######Make panel
#####################################################################################################################
#####put them all in the same pdf
#really helpful: http://www.sthda.com/english/rpkgs/ggpubr/index.html
#install.packages("ggpubr")
library(ggpubr)

#make one with a legend
alphadiv <- ggarrange(bacscatter, fungalscatter, ncol = 2, common.legend=TRUE,legend="bottom")
alphadiv

#make one without a legend
alphadiv_nolegend <- ggarrange(bacscatter2, fungalscatter2)
alphadiv_nolegend 

pdf("Figures/zotuvotu_richness/Figure1_alpharichness_updated.pdf", height=6, width=10)
alphadiv
dev.off()

################################################################################################
######Fungi: avg.dist bray-curtis dissimilarity matrices - median and square root transformed for T1 T2 T3 trasplant
################################################################################################
#bray curtis
fungi_T1_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T1_otu.csv", row.names=1 )
fungi_T1_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T1_zotu.csv", row.names=1 )

fungi_T2_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T2_otu.csv", row.names=1 )
fungi_T2_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T2_zotu.csv", row.names=1 )

fungi_T3_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T3_otu.csv", row.names=1 )
fungi_T3_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T3_zotu.csv", row.names=1 )

################################################################################################
######Bacteria: avg.dist bray-curtis dissimilarity matrices - median and square root transformed for T1 T2 T3 trasplant
################################################################################################
#bray-curtis
bac_T1_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T1_otu.csv", row.names=1 )
bac_T1_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T1_zotu.csv", row.names=1 )

bac_T2_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T2_otu.csv", row.names=1 )
bac_T2_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T2_zotu.csv", row.names=1 )

bac_T3_bray_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T3_otu.csv", row.names=1 )
bac_T3_bray_zotu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T3_zotu.csv", row.names=1 )


################################################################################################
######Make everything vectors and dataframes that can be put into ggplot2 figure
################################################################################################
#T1
#make this a vector so I can make a ggplot out of it
fungi_T1_ESV_1 <- as.dist(fungi_T1_bray_zotu)
fungi_T1_ESV <- fungi_T1_ESV_1 [lower.tri(fungi_T1_ESV_1)]

fungi_T1_OTU_1 <- as.dist(fungi_T1_bray_otu)
fungi_T1_OTU<- fungi_T1_OTU_1 [lower.tri(fungi_T1_OTU_1)]

#make it a dataframe
fungiT1<- as.data.frame(cbind(fungi_T1_ESV,fungi_T1_OTU))
fungiT1$Timepoint <- rep("T1",nrow(fungiT1))

#T2
#make this a vector so I can make a ggplot out of it
fungi_T2_ESV_1 <- as.dist(fungi_T2_bray_zotu)
fungi_T2_ESV <-fungi_T2_ESV_1[lower.tri(fungi_T2_ESV_1)]

fungi_T2_OTU_1 <- as.dist(fungi_T2_bray_otu)
fungi_T2_OTU<- fungi_T2_OTU_1 [lower.tri(fungi_T2_OTU_1)]

#make it a dataframe
fungiT2<- as.data.frame(cbind(fungi_T2_ESV,fungi_T2_OTU))
fungiT2$Timepoint <- rep("T2",nrow(fungiT2))

#T3
#make this a vector so I can make a ggplot out of it
fungi_T3_ESV_1 <- as.dist(fungi_T3_bray_zotu)
fungi_T3_ESV <- fungi_T3_ESV_1[lower.tri(fungi_T3_ESV_1)]

fungi_T3_OTU_1 <- as.dist(fungi_T3_bray_otu)
fungi_T3_OTU<- fungi_T3_OTU_1[lower.tri(fungi_T3_OTU_1)]

#make it a dataframe
fungiT3<- as.data.frame(cbind(fungi_T3_ESV,fungi_T3_OTU))
fungiT3$Timepoint <- rep("T3",nrow(fungiT3))

#combine all
names(fungiT1) <- c("ESV","OTU","Timepoint")
names(fungiT2)<- c("ESV","OTU","Timepoint")
names(fungiT3)<- c("ESV","OTU","Timepoint")
fungi_T1T2T3_beta <- rbind(fungiT1, fungiT2, fungiT3)

fungi_T1T2T3_beta$Timepoint <- as.factor(fungi_T1T2T3_beta$Timepoint)
################################################################################################
######bacteria Make everything vectors and dataframes that can be put into ggplot2 figure
################################################################################################
#T1
#make this a vector so I can make a ggplot out of it
bac_T1_ESV_1 <- as.dist(bac_T1_bray_zotu)
bac_T1_ESV <- bac_T1_ESV_1[lower.tri(bac_T1_ESV_1)]

bac_T1_OTU_1 <- as.dist(bac_T1_bray_otu)
bac_T1_OTU<- bac_T1_OTU_1 [lower.tri(bac_T1_OTU_1)]

#make it a dataframe
bacT1<- as.data.frame(cbind(bac_T1_ESV,bac_T1_OTU))
bacT1$Timepoint <- rep("T1",nrow(bacT1))

#T2
#make this a vector so I can make a ggplot out of it
bac_T2_ESV_1 <- as.dist(bac_T2_bray_zotu)
bac_T2_ESV <- bac_T2_ESV_1 [lower.tri(bac_T2_ESV_1)]

bac_T2_OTU_1 <- as.dist(bac_T2_bray_otu)
bac_T2_OTU<- bac_T2_OTU_1 [lower.tri(bac_T2_OTU_1)]

#make it a dataframe
bacT2<- as.data.frame(cbind(bac_T2_ESV,bac_T2_OTU))
bacT2$Timepoint <- rep("T2",nrow(bacT2))

#T3
#make this a vector so I can make a ggplot out of it
bac_T3_ESV_1 <- as.dist(bac_T3_bray_zotu)
bac_T3_ESV <-bac_T3_ESV_1 [lower.tri(bac_T3_ESV_1)]

bac_T3_OTU_1 <- as.dist(bac_T3_bray_otu)
bac_T3_OTU<- bac_T3_OTU_1[lower.tri(bac_T3_OTU_1)]

#make it a dataframe
bacT3<- as.data.frame(cbind(bac_T3_ESV,bac_T3_OTU))
bacT3$Timepoint <- rep("T3",nrow(bacT3))

#combine all
names(bacT1) <- c("ESV","OTU","Timepoint")
names(bacT2)<- c("ESV","OTU","Timepoint")
names(bacT3)<- c("ESV","OTU","Timepoint")
bac_T1T2T3_beta <- rbind(bacT1, bacT2, bacT3)

bac_T1T2T3_beta$Timepoint <- as.factor(bac_T1T2T3_beta$Timepoint)
################################################################################################
######Make a  beta-diversity figure plotting dissimilarity matrices against each other in ggplot2
################################################################################################
#change the levels from T1 T2 T3 to 6,12,18 months
levels(bac_T1T2T3_beta$Timepoint) <- c("6 months","12 months","18 months")

bacterialscatter_beta <- ggscatter(bac_T1T2T3_beta, x = "OTU", y = "ESV",
                                   add = "reg.line",                         # Add regression line
                                   conf.int = TRUE,                          # Add confidence interval
                                   color = "Timepoint",            # Color by groups "cyl"
                                   shape = "Timepoint",                             # Change point shape by groups "cyl"
                                   size=1,                          # Change size to 2
                                   alpha=0.5                        # Add transparency to points
)+
  
  #stat_cor(aes(color = Timepoint),label.x = 100,label.y.npc=0.1, method="pearson")   + # Add correlation coefficient
  labs(x="Bacterial OTU beta-diversity",y="Bacterial ESV beta-diversity", title="C") + # change labels
  theme(axis.title=element_text(size=14), axis.text.x=element_text(size=14, angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger  


bacterialscatter_beta 

################################################################################################
######Make a  beta-diversity figure plotting dissimilarity matrices against each other in ggplot2
################################################################################################
#change the levels from T1 T2 T3 to 6,12,18 months
levels(fungi_T1T2T3_beta$Timepoint) <- c("6 months","12 months","18 months")

fungalscatter_beta <- ggscatter(fungi_T1T2T3_beta, x = "OTU", y = "ESV",
                           add = "reg.line",                         # Add regression line
                           conf.int = TRUE,                          # Add confidence interval
                           color = "Timepoint",            # Color by groups "cyl"
                           shape = "Timepoint",                             # Change point shape by groups "cyl"
                           size=1 ,                         # Change size to 2
                           alpha=0.5                        # Add transparency to points
)+
  #stat_cor(aes(color = Timepoint),label.x = 100,label.y.npc=0.1, method="pearson")   + # Add correlation coefficient
  labs(x="Fungal OTU beta-diversity",y="Fungal ESV beta-diversity", title="D") + # change labels
  theme(axis.title=element_text(size=14), axis.text.x=element_text(size=14, angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger  
        
fungalscatter_beta 


################################################################################################
######Make a  beta-diversity figure plotting dissimilarity matrices against each other
################################################################################################
betadiv <- ggarrange(bacterialscatter_beta, fungalscatter_beta, ncol = 2, common.legend=TRUE,legend="bottom")


#make versions without legends
fungalscatter_beta2 <- fungalscatter_beta  +   theme(legend.position="none")#remove legend
bacterialscatter_beta2  <- bacterialscatter_beta  +   theme(legend.position="none")#remove legend
#arrange into one panel
betadiv_nolegend <-  ggarrange(bacterialscatter_beta2,fungalscatter_beta2)
betadiv_nolegend                              
                              
#arrange alpha nad beta to one
alpha_beta_div <- ggarrange(alphadiv, betadiv_nolegend , nrow=2, common.legend=TRUE,legend="bottom")
alpha_beta_div

pdf("Figures/zotuvotu_richness/Figure1_alphabeta.pdf", height=10, width=10)
alpha_beta_div
dev.off()

