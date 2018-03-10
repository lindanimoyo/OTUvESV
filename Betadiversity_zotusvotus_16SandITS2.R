#December 12, 2017
#upload all OTU tables and all ZOTU tables
#make bray curtis matrices
#do mantel loop to get an output of all of the R values and p values
#make a dissimilarity matrix mantel correlation figure for a fungi and bacteria

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

#source in functions
source('~/Dropbox/StatsandProgramming/source/gettaxondd.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/getrowsums.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/pvalueadonis.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/pvaluelegend.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/avgdist2.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/getbraydistfunction.R', chdir = TRUE)


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

#jaccard
fungi_T1_jac_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T1_otu_jac.csv",row.names=1 )
fungi_T1_jac_zotu <-  read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T1_zotu_jac.csv",,row.names=1 )

fungi_T2_jac_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T2_otu_jac.csv",row.names=1 )
fungi_T2_jac_zotu <-  read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T2_zotu_jac.csv",,row.names=1 )

fungi_T3_jac_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T3_otu_jac.csv",row.names=1 )
fungi_T3_jac_zotu <-  read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_ITS2_T3_zotu_jac.csv",,row.names=1 )

#check that row and column names all match
row.names(fungi_T1_bray_otu) == row.names(fungi_T1_bray_zotu)
row.names(fungi_T1_bray_otu) ==colnames(fungi_T1_bray_otu)
colnames(fungi_T1_bray_otu) ==colnames(fungi_T1_bray_zotu)

row.names(fungi_T2_bray_otu) == row.names(fungi_T2_bray_zotu)
row.names(fungi_T2_bray_otu) ==colnames(fungi_T2_bray_otu)
colnames(fungi_T2_bray_otu) ==colnames(fungi_T2_bray_zotu)

row.names(fungi_T3_bray_otu) == row.names(fungi_T3_bray_zotu)
row.names(fungi_T3_bray_otu) ==colnames(fungi_T3_bray_otu)
colnames(fungi_T3_bray_otu) ==colnames(fungi_T3_bray_zotu)

row.names(fungi_T1_jac_otu) == row.names(fungi_T1_jac_zotu)
row.names(fungi_T1_jac_otu) ==colnames(fungi_T1_jac_otu)
colnames(fungi_T1_jac_otu) ==colnames(fungi_T1_jac_zotu)

row.names(fungi_T2_jac_otu) == row.names(fungi_T2_jac_zotu)
row.names(fungi_T2_jac_otu) ==colnames(fungi_T2_jac_otu)
colnames(fungi_T2_jac_otu) ==colnames(fungi_T2_jac_zotu)

row.names(fungi_T3_jac_otu) == row.names(fungi_T3_jac_zotu)
row.names(fungi_T3_jac_otu) ==colnames(fungi_T3_jac_otu)
colnames(fungi_T3_jac_otu) ==colnames(fungi_T3_jac_zotu)
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

#jaccard
bac_T1_jac_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T1_otu_jac.csv",row.names=1 )
bac_T1_jac_zotu <-  read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T1_zotu_jac.csv",row.names=1 )

bac_T2_jac_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T2_otu_jac.csv",row.names=1 )
bac_T2_jac_zotu <-  read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T2_zotu_jac.csv",,row.names=1 )

bac_T3_jac_otu <- read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T3_otu_jac.csv",row.names=1 )
bac_T3_jac_zotu <-  read.csv("data/DissMatricesforPRIMER/braycurtis_rarefied_dm_sqrt_16S_T3_zotu_jac.csv",,row.names=1 )

#check that row and column names all match
row.names(bac_T1_bray_otu) == row.names(bac_T1_bray_zotu)
row.names(bac_T1_bray_otu) ==colnames(bac_T1_bray_otu)
colnames(bac_T1_bray_otu) ==colnames(bac_T1_bray_zotu)

row.names(bac_T2_bray_otu) == row.names(bac_T2_bray_zotu)
row.names(bac_T2_bray_otu) ==colnames(bac_T2_bray_otu)
colnames(bac_T2_bray_otu) ==colnames(bac_T2_bray_zotu)

row.names(bac_T3_bray_otu) == row.names(bac_T3_bray_zotu)
row.names(bac_T3_bray_otu) ==colnames(bac_T3_bray_otu)
colnames(bac_T3_bray_otu) ==colnames(bac_T3_bray_zotu)

row.names(bac_T1_jac_otu) == row.names(bac_T1_jac_zotu)
row.names(bac_T1_jac_otu) ==colnames(bac_T1_jac_otu)
colnames(bac_T1_jac_otu) ==colnames(bac_T1_jac_zotu)

row.names(bac_T2_jac_otu) == row.names(bac_T2_jac_zotu)
row.names(bac_T2_jac_otu) ==colnames(bac_T2_jac_otu)
colnames(bac_T2_jac_otu) ==colnames(bac_T2_jac_zotu)

row.names(bac_T3_jac_otu) == row.names(bac_T3_jac_zotu)
row.names(bac_T3_jac_otu) ==colnames(bac_T3_jac_otu)
colnames(bac_T3_jac_otu) ==colnames(bac_T3_jac_zotu)

################################################################################################
######Run Mantel tests to Mantel r and pvalues
################################################################################################

manteltest_bac_T1 <- mantel(bac_T1_bray_otu,bac_T1_bray_zotu, permutations=999)
manteltest_bac_T2 <- mantel(bac_T2_bray_otu,bac_T2_bray_zotu, permutations=999)
manteltest_bac_T3 <- mantel(bac_T3_bray_otu,bac_T3_bray_zotu, permutations=999)

manteltest_bac_T1_jac <- mantel(bac_T1_jac_otu,bac_T1_jac_zotu, permutations=999)
manteltest_bac_T2_jac <- mantel(bac_T2_jac_otu,bac_T2_jac_zotu, permutations=999)
manteltest_bac_T3_jac <- mantel(bac_T3_jac_otu,bac_T3_jac_zotu, permutations=999)

manteltest_fungi_T1 <- mantel(fungi_T1_bray_otu,fungi_T1_bray_zotu, permutations=999)
manteltest_fungi_T2 <- mantel(fungi_T2_bray_otu,fungi_T2_bray_zotu, permutations=999)
manteltest_fungi_T3 <- mantel(fungi_T3_bray_otu,fungi_T3_bray_zotu, permutations=999)

manteltest_fungi_T1_jac <- mantel(fungi_T1_jac_otu,fungi_T1_jac_zotu, permutations=999)
manteltest_fungi_T2_jac <- mantel(fungi_T2_jac_otu,fungi_T2_jac_zotu, permutations=999)
manteltest_fungi_T3_jac <- mantel(fungi_T3_jac_otu,fungi_T3_jac_zotu, permutations=999)

#make dataframe of mantel r value results for otu v zotu 
mantelresults <- as.data.frame(rbind(cbind("Bacteria","T1",manteltest_bac_T1$statistic,manteltest_bac_T1$signif,manteltest_bac_T1_jac$statistic,manteltest_bac_T1_jac$signif),
                                    cbind("Bacteria","T2",manteltest_bac_T2$statistic,manteltest_bac_T2$signif,manteltest_bac_T2_jac$statistic,manteltest_bac_T2_jac$signif),
                                    cbind("Bacteria","T3",manteltest_bac_T3$statistic,manteltest_bac_T3$signif,manteltest_bac_T3_jac$statistic,manteltest_bac_T3_jac$signif),
                                    cbind("Fungi","T1",manteltest_fungi_T1$statistic,manteltest_fungi_T1$signif,manteltest_fungi_T1_jac$statistic,manteltest_fungi_T1_jac$signif),
                                    cbind("Fungi","T2",manteltest_fungi_T2$statistic,manteltest_fungi_T2$signif,manteltest_fungi_T2_jac$statistic,manteltest_fungi_T2_jac$signif),
                                    cbind("Fungi","T3",manteltest_fungi_T3$statistic,manteltest_fungi_T3$signif,manteltest_fungi_T3_jac$statistic,manteltest_fungi_T3_jac$signif)))
#add names to dataframe
names(mantelresults)<- c("Microbe","Timepoint","Bray_Curtis_Mantel_r","Bray_Curtis_P-val","Jaccard_Mantel_r","Jaccard_P-val")

#make mantel results numeric and only 2 decimals
mantelresults$Bray_Curtis_Mantel_r <- as.numeric(as.character(mantelresults$Bray_Curtis_Mantel_r))
mantelresults$Bray_Curtis_Mantel_r <- round(mantelresults$Bray_Curtis_Mantel_r, digits=2)
mantelresults$Jaccard_Mantel_r <- as.numeric(as.character(mantelresults$Jaccard_Mantel_r))
mantelresults$Jaccard_Mantel_r <- round(mantelresults$Jaccard_Mantel_r, digits=2)
#show data frame
mantelresults

write.csv(mantelresults, "results/betadiversity_otuvzotu_mantelresults.csv")
################################################################################################
######getting equations of the lines
################################################################################################
fit_bac_T1  <- lm(as.dist(bac_T1_bray_zotu) ~ as.dist(bac_T1_bray_otu))
fit_bac_T2  <- lm(as.dist(bac_T2_bray_zotu) ~ as.dist(bac_T2_bray_otu))
fit_bac_T3  <- lm(as.dist(bac_T3_bray_zotu) ~ as.dist(bac_T3_bray_otu))

fit_bac_T1_jac  <- lm(as.dist(bac_T1_jac_zotu) ~ as.dist(bac_T1_jac_otu))
fit_bac_T2_jac  <- lm(as.dist(bac_T2_jac_zotu) ~ as.dist(bac_T2_jac_otu))
fit_bac_T3_jac  <- lm(as.dist(bac_T3_jac_zotu) ~ as.dist(bac_T3_jac_otu))

fit_fungi_T1  <- lm(as.dist(fungi_T1_bray_zotu) ~ as.dist(fungi_T1_bray_otu))
fit_fungi_T2  <- lm(as.dist(fungi_T2_bray_zotu) ~ as.dist(fungi_T2_bray_otu))
fit_fungi_T3  <- lm(as.dist(fungi_T3_bray_zotu) ~ as.dist(fungi_T3_bray_otu))

fit_fungi_T1_jac  <- lm(as.dist(fungi_T1_jac_zotu) ~ as.dist(fungi_T1_jac_otu))
fit_fungi_T2_jac  <- lm(as.dist(fungi_T2_jac_zotu) ~ as.dist(fungi_T2_jac_otu))
fit_fungi_T3_jac  <- lm(as.dist(fungi_T3_jac_zotu) ~ as.dist(fungi_T3_jac_otu))

summary(fit_bac_T1)
summary(fit_bac_T2)
summary(fit_bac_T3)

summary(fit_bac_T1_jac)
summary(fit_bac_T2_jac)
summary(fit_bac_T3_jac)

summary(fit_fungi_T1)
summary(fit_fungi_T2)
summary(fit_fungi_T3)

summary(fit_fungi_T1_jac)
summary(fit_fungi_T2_jac)
summary(fit_fungi_T3_jac)
###########bray-curtis
equationofline_results <- rbind(cbind("Bacteria","T1",coef(fit_bac_T1)[1],coef(fit_bac_T1)[2],summary(fit_bac_T1)$r.squared),
                                     cbind("Bacteria","T2",coef(fit_bac_T2)[1],coef(fit_bac_T2)[2],summary(fit_bac_T2)$r.squared),
                                     cbind("Bacteria","T3",coef(fit_bac_T3)[1],coef(fit_bac_T3)[2],summary(fit_bac_T3)$r.squared),
                                     cbind("Fungi","T1",coef(fit_fungi_T1)[1],coef(fit_fungi_T1)[2],summary(fit_fungi_T1)$r.squared),
                                     cbind("Fungi","T2",coef(fit_fungi_T2)[1],coef(fit_fungi_T2)[2],summary(fit_fungi_T2)$r.squared),
                                     cbind("Fungi","T3",coef(fit_fungi_T3)[1],coef(fit_fungi_T3)[2],summary(fit_fungi_T3)$r.squared))

row.names(equationofline_results ) <- c(1,2,3,4,5,6)
equationofline_results <- as.data.frame(equationofline_results)
names(equationofline_results)<- c("Microbe","Timeopoint","Intercept","Slope","R squared")
equationofline_results$Intercept <- as.numeric(as.character(equationofline_results$Intercept))
equationofline_results$Intercept <- round(equationofline_results$Intercept, digits=2)
equationofline_results$Slope <- as.numeric(as.character(equationofline_results$Slope))
equationofline_results$Slope <- round(equationofline_results$Slope, digits=2)
equationofline_results$"R squared" <- as.numeric(as.character(equationofline_results$"R squared"))
equationofline_results$"R squared"<- round(equationofline_results$"R squared", digits=2)

equationofline_results$Metric <- rep("Bray-Curtis", nrow(equationofline_results))
equationofline_results #make this for alpha diversity too

###########jaccard
equationofline_results_jac <- rbind(cbind("Bacteria","T1",coef(fit_bac_T1_jac)[1],coef(fit_bac_T1_jac)[2],summary(fit_bac_T1_jac)$r.squared),
                                cbind("Bacteria","T2",coef(fit_bac_T2_jac)[1],coef(fit_bac_T2_jac)[2],summary(fit_bac_T2_jac)$r.squared),
                                cbind("Bacteria","T3",coef(fit_bac_T3_jac)[1],coef(fit_bac_T3_jac)[2],summary(fit_bac_T3_jac)$r.squared),
                                cbind("Fungi","T1",coef(fit_fungi_T1_jac)[1],coef(fit_fungi_T1_jac)[2],summary(fit_fungi_T1_jac)$r.squared),
                                cbind("Fungi","T2",coef(fit_fungi_T2_jac)[1],coef(fit_fungi_T2_jac)[2],summary(fit_fungi_T2_jac)$r.squared),
                                cbind("Fungi","T3",coef(fit_fungi_T3_jac)[1],coef(fit_fungi_T3_jac)[2],summary(fit_fungi_T3_jac)$r.squared))

row.names(equationofline_results_jac ) <- c(1,2,3,4,5,6)
equationofline_results_jac <- as.data.frame(equationofline_results_jac)
names(equationofline_results_jac)<- c("Microbe","Timeopoint","Intercept","Slope","R squared")
equationofline_results_jac$Intercept <- as.numeric(as.character(equationofline_results$Intercept))
equationofline_results_jac$Intercept <- round(equationofline_results_jac$Intercept, digits=2)
equationofline_results_jac$Slope <- as.numeric(as.character(equationofline_results_jac$Slope))
equationofline_results_jac$Slope <- round(equationofline_results_jac$Slope, digits=2)
equationofline_results_jac$"R squared" <- as.numeric(as.character(equationofline_results_jac$"R squared"))
equationofline_results_jac$"R squared"<- round(equationofline_results_jac$"R squared", digits=2)

equationofline_results_jac$Metric <- rep("Jaccard", nrow(equationofline_results_jac))
equationofline_results_jac #make this for alpha diversity too

betadiversityresults <- rbind(equationofline_results, equationofline_results_jac)

betadiversityresults
write.csv(betadiversityresults, "results/betadiversitycomparisons_otuvzotu_fungiandbacteria.csv")
#also do with jaccard for comparison
################################################################################################
######Make a figure plotting dissimilarity matrices against each other
################################################################################################
shapesfortimepoint <- c(16, 17, 15)
library(scales)
pdf("Figures/otuvzotu_dissimilarity/BetadiversityFigure_otuvzotu_fungiandbac.pdf", height=6, width=10)
par(mfrow=c(1,2))
##Bacteria
 #mantel test for T1 bacteria: make ZOTU (ESVS) and on y axis
manteltest <- mantel(x=bac_T1_bray_otu,y=bac_T1_bray_zotu, permutations=999) #run Mantel test
plot(y=as.dist(bac_T1_bray_zotu),x=as.dist(bac_T1_bray_otu), #plot zotu on x axis, otu on y axis
     cex.lab=1, #make labels larger
     pch=16, #shape from ggplot
     col=alpha("#F8766D",0.8), #color from ggplot
     xlab="Bacterial OTU Beta-diversity", ylab="Bacterial ESV Beta-Diversity", #fix x and y labels and do ESV instead ot OTU
     ylim=c(0,1),xlim=c(0,1)) #make x and y lim the same
#pvaluelegend(manteltest$statistic, manteltest$signif) #show legend
abline(lm(as.dist(bac_T1_bray_zotu) ~as.dist(bac_T1_bray_otu)),col="#F8766D",lwd=2) #plot equation of the line, make same colors as ggplot and thicker line width like ggplot

#add test for T2
par(new=TRUE) #make it so next plot adds on to current plot
plot(as.dist(bac_T2_bray_zotu) ~ as.dist(bac_T2_bray_otu), col=alpha("#009E73",0.5),ylim=c(0,1),xlim=c(0,1),axes=FALSE,  xlab="", ylab="", pch=17,cex.lab=1 )
abline(lm(as.dist(bac_T2_bray_zotu) ~as.dist(bac_T2_bray_otu)),col="#009E73",lwd=2)

#add test for T3
par(new=TRUE) #make it so next plot adds on to current plot
plot(as.dist(bac_T3_bray_zotu) ~ as.dist(bac_T3_bray_otu), col=alpha("#56B4E9",0.1), ylim=c(0,1),xlim=c(0,1),axes=FALSE,  xlab="", ylab="", pch=15,cex.lab=1)
abline(lm(as.dist(bac_T3_bray_zotu) ~as.dist(bac_T3_bray_otu)),col="#56B4E9",lwd=2)

mtext("A",line=1, adj=0, cex=1.5)

##FUngi
#mantel test for T1 fungi 
manteltest1 <- mantel(fungi_T1_bray_otu,fungi_T1_bray_zotu, permutations=999) 
#plot zotu as y against otu as x, same colors as ggplot2, same axes for all, change ZOTU to ESV in text
plot(as.dist(fungi_T1_bray_zotu)~as.dist(fungi_T1_bray_otu), 
     cex.lab=1,
     col=alpha("#F8766D",0.8),
     pch=16, #shape from ggplot
     xlab="Fungal OTU Beta-diversity", ylab="Fungal ESV Beta-diversity",ylim=c(0,1),xlim=c(0,1))
abline(lm(as.dist(fungi_T1_bray_zotu) ~as.dist(fungi_T1_bray_otu)),col="#F8766D",lwd=2)
#pvaluelegend(manteltest1$statistic, manteltest1$signif)
#add test for T2
par(new=TRUE) #make it so next plot adds on to current plot
plot(as.dist(fungi_T2_bray_zotu) ~ as.dist(fungi_T2_bray_otu), col=alpha("#009E73",0.5),ylim=c(0,1),xlim=c(0,1),axes=FALSE,  xlab="", ylab="", pch=17,cex.lab=1)
abline(lm(as.dist(fungi_T2_bray_zotu) ~as.dist(fungi_T2_bray_otu)),col="#009E73",lwd=2)

#add test for T3
par(new=TRUE) #make it so next plot adds on to current plot
plot(as.dist(fungi_T3_bray_zotu) ~ as.dist(fungi_T3_bray_otu), 
     col=alpha("#56B4E9",0.1), #make more transparent 
     ylim=c(0,1),xlim=c(0,1),axes=FALSE,  xlab="", ylab="", pch=15, cex=1)
abline(lm(as.dist(fungi_T3_bray_zotu) ~as.dist(fungi_T3_bray_otu)),col="#56B4E9",lwd=2)
mtext("B",line=1, adj=0, cex=1.5)

#same colors as ggplot alpha diversity plot
colors <- c("#F8766D","#009E73","#56B4E9")
legend("bottomright", legend=c("6 months","12 months","18 months"), col=colors, bty="n", pch=shapesfortimepoint, cex=1.5)

dev.off()


######################################################################################################################
###### Make example NMDS plot and stress plots - Bac and Fungi, OTU and ZOTU
######################################################################################################################
#make the nmds 

plotnmds_fungi_T3_otu <- metaMDS(fungi_T3_bray_otu, k=2, trymax=100)
plotnmds_fungi_T3_otu 
stressplot(plotnmds_fungi_T3_otu)

plotnmds_fungi_T3_zotu <- metaMDS(fungi_T3_bray_zotu, k=2, trymax=100)
plotnmds_fungi_T3_zotu
stressplot(plotnmds_fungi_T3_zotu)

plotnmds_bac_T3_otu <- metaMDS(bac_T3_bray_otu, k=3, trymax=100)
plotnmds_bac_T3_otu
stressplot(plotnmds_bac_T3_otu)

plotnmds_bac_T3_zotu <- metaMDS(bac_T3_bray_zotu, k=3, trymax=100)
plotnmds_bac_T3_zotu
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

pdf("Figures/otuvzotu_dissimilarity/NMDS_ITS2_T3_forpub.pdf",height=6, width=10)
par(mfrow=c(1,2), mai=c(1,0.8,1,0.2))
plot(scores(plotnmds_fungi_T3_otu), col=map_T3_fungi$colorbyinoc, pch=map_T3_fungi$shapesbysite, cex=1, ylim=c(-0.4,0.6), xlim=c(-0.5,0.6))
mtext("A Fungal OTU",line=1, adj=0, cex=2)
plot(scores(plotnmds_fungi_T3_zotu), col=map_T3_fungi$colorbyinoc, pch=map_T3_fungi$shapesbysite, cex=1,ylim=c(-0.4,0.6), xlim=c(-0.5,0.6))
mtext("B Fungal ESV",line=1, adj=0, cex=2)
dev.off()

pdf("Figures/otuvzotu_dissimilarity/NMDS_16S_T3_forpub.pdf",height=6, width=10)
par(mfrow=c(1,2), mai=c(1,0.8,1,0.2))
plot(scores(plotnmds_bac_T3_otu), col=map_T3_bac$colorbysite, pch=map_T3_bac$ShapesbyInoculum, cex=1, ylim=c(-0.4,0.4), xlim=c(-0.55,0.5))
mtext("A Bacterial OTU",line=1, adj=0, cex=2)
plot(scores(plotnmds_bac_T3_zotu), col=map_T3_bac$colorbysite, pch=map_T3_bac$ShapesbyInoculum, cex=1,ylim=c(-0.4,0.4),xlim=c(-0.55,0.5))
mtext("B Bacterial ESV",line=1, adj=0, cex=2)
dev.off()

par(mfrow=c(1,1),par(mai=c(1.02,0.82,0.82,0.42)))

###########################################
###### Put all four on the same figure 
############################################

#make vectors for legend
sitenames <- c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine")
sitecolors  <- c("red","orange","green","blue","purple")
shapes <- c(17,3,16,15,18)

#mai  c(bottom, left, top, right)
#side	
#on which side of the plot (1=bottom, 2=left, 3=top, 4=right).
plot.new()
pdf("Figures/otuvzotu_dissimilarity/NMDS_ITS2and16S_T3_otuvzotu.pdf", height=8, width=8)
par(mfrow=c(2,2)) #make 2 columns by 2 rows
#Bacteria 
#Panel A
par(mar=c(0.2,5,5,0.2))
plot(scores(plotnmds_bac_T3_otu),col= map_T3_bac$colorbysite, pch=16, cex =1.5, ylim=c(-0.55, 0.55),  xlim=c(-0.55, 0.55), xaxt='n')#remove x axis
mtext("OTU", line=1, adj=0.5, cex=2, side=3)
mtext("A", line=-2, adj=0.05, cex=2, side=3)
#mtext("16S; Stress = 0.067", line=-2, adj=0.9, cex=1, side=3)
#plotnmds_bac_T3_otu$stress
mtext("Stress = 0.067", line=-2, adj=0.9, cex=1, side=3)

#Panel B
par(mar=c(0.2,0.2,5,5))
plot(scores(plotnmds_bac_T3_zotu),col= map_T3_bac$colorbysite, pch=16, cex =1.5, ylim=c(-0.55, 0.55),  xlim=c(-0.55, 0.55),xaxt='n', yaxt='n', ylab="") #remove x and y axis
mtext("ESV", line=1, adj=0.5, cex=2, side=3)
mtext("B", line=-2, adj=0.05, cex=2,side=3)
#plotnmds_bac_T3_zotu$stress
mtext("Stress = 0.078", line=-2, adj=0.9, cex=1, side=3)
mtext("Bacteria", line=2, adj=0.5, cex=2, side=4)

#Fungi
#Panel C
par(mar=c(5,5,0.2,0.2))
plot(scores(plotnmds_fungi_T3_otu),col= map_T3_fungi$colorbyinoc, pch=16, cex =1.5,ylim=c(-0.55, 0.55), xlim=c(-0.55, 0.55)) 
mtext("C", line=-2, adj=0.05, cex=2,side=3)
#plotnmds_fungi_T3_otu$stress
mtext("Stress = 0.15", line=-2, adj=0.9, cex=1, side=3)
#legend("bottomleft",legend=sitenames,col=sitecolors, pch=16,bty="l",cex=0.8)

#Panel D
par(mar=c(5,0.2,0.2,5))
plot(scores(plotnmds_fungi_T3_zotu),col= map_T3_fungi$colorbyinoc, pch=16, cex =1.5, ylim=c(-0.55, 0.55), xlim=c(-0.55, 0.55),yaxt='n', ylab="") #remove y axis
mtext("D", line=-2, adj=0.05, cex=2,side=3)
#plotnmds_fungi_T3_zotu$stress
mtext("Stress = 0.14", line=-2, adj=0.9, cex=1, side=3)
mtext("Fungi", line=2, adj=0.5, cex=2, side=4)
#legend("bottomright",legend=sitenames,pch=shapes,col="black", bty="l", cex=0.8) #one legend for whole thing, combine colors and shapes

dev.off()


################################################################################################
######try to make it a loop but this doesnt work
################################################################################################
##this loop ins't working for some reason
#create lists of dissimilarity matrices for otu v zotu
mylist_otu <-list(bac_T1_bray_otu,bac_T2_bray_otu,bac_T3_bray_otu,fungi_T1_bray_otu,fungi_T2_bray_otu,fungi_T3_bray_otu)
mylist_zotu <-list(bac_T1_bray_zotu,bac_T2_bray_zotu,bac_T3_bray_zotu,fungi_T1_bray_zotu,fungi_T2_bray_zotu,fungi_T3_bray_zotu)

#create empty vectors to store the data
mantelresults <- NULL
P.val <- NULL
R2 <- NULL


#run for loop to do mantel correlations
for(k in 1:length(mylist_otu)){
  for(j in 1:length(mylist_zotu)){
    #do the adonis test on the beta diversity against tree
    manteltest <- mantel(mylist_otu[[k]]~mylist_zotu[[j]], permutations=999)
    #create a vector of values of the R2 values of the tree host effect
    R2 <- rbind(R2, manteltest$statistic)
    #create a vector of values of the p-values of the tree hose effect
    P.val <- rbind(P.val,manteltest$signif)
  }
}
#create a matrix of the R2 values
mantelresults  <-  rbind(mantelresults ,R2)
#add columns of the pvalues
mantelresults <- cbind(mantelresults ,P.val)
#check out the matrix
mantelresults 
#to create a data frame with the names of the distance matrices used and the results, create a vector
#of the names of distance matrices and rbind it to adonisresults
mylist_names<-c("bac_T1_bray","bac_T2_bray","bac_T3_bray","fungi_T1_bray",
                "fungi_T2_bray","fungi_T3_bray")
#bind the names to my results matrix
mantelresults_names <- cbind(mylist_names,mantelresults )
#now look at my results with names
mantelresults 





