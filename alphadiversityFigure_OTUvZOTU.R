#Dec 11, 2017

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
correlations <- read.csv("data/correlationsforotuvzotu.csv", header=TRUE)
fungalrichness <- read.csv("data/ITS2_OTUvZOTUall.csv", row.names=1)

################################################################################################
######getting equations of the lines
################################################################################################
#subset bacteria
bacrichness_T1 <- bacrichness[bacrichness$Timepoint=="T1", ]
bacrichness_T2 <- bacrichness[bacrichness$Timepoint=="T2", ]
bacrichness_T3 <- bacrichness[bacrichness$Timepoint=="T3", ]

#subset fungi
fungirichness_T1 <- fungalrichness [fungalrichness $Timepoint=="T1", ]
fungirichness_T2 <- fungalrichness [fungalrichness $Timepoint=="T2", ]
fungirichness_T3 <- fungalrichness [fungalrichness $Timepoint=="T3", ]

########################################################################################################################
######Make loop to get equation of line results intercept, slope, R2 and p value for linear fit of OTU to Zotu richness
#########################################################################################################################

#create a list of richness data frames for bacteria and fungi for each time point 
mylist_otu <-list(bacrichness_T1,bacrichness_T2,bacrichness_T3,fungirichness_T1,fungirichness_T2,fungirichness_T3)
bacrichness_T1[9]
bacrichness_T1[24]
head(fungirichness_T1[9])
head(fungirichness_T1[24])
head(bacrichness_T1[11])
head(bacrichness_T1[26])
head(fungirichness_T1[9])
head(fungirichness_T1[24])
#names(bacrichness_T1)     
#richness <- 9 vs 24
#chao1 <- 3 vs 18
#simpson <- 11 v 26
#shannon <- 12 v 27
#berger_parker <- 1, 16


#create empty vectors to store the data
lmresults <- NULL
Intercept <- NULL
Slope <- NULL
P.val <- NULL
R2 <- NULL
P.val <- NULL

#run for loop to do linear models
for(k in 1:length(mylist_otu)){
    #do the linear model test of otu richness against zotu richness
    fit_richness <- lm(mylist_otu[[k]][[9]]~mylist_otu[[k]][[24]])
    #do the linear model test of otu richness against zotu chao1
    fit_chao1<- lm(mylist_otu[[k]][[3]]~mylist_otu[[k]][[18]])
    #do the linear model test of otu richness against zotu Simpson
    fit_simpson <- lm(mylist_otu[[k]][[11]]~mylist_otu[[k]][[26]])
    #do the linear model test of otu richness against zotu Shannon
    fit_shannon <- lm(mylist_otu[[k]][[12]]~mylist_otu[[k]][[27]])
    #do the linear model test of otu richness against zotu berger parker
    fit_bp <- lm(mylist_otu[[k]][[1]]~mylist_otu[[k]][[16]])
    #create a vector of values of the intercept
    Intercept <- rbind(Intercept, coef(fit_richness)[1],coef(fit_chao1)[1],coef(fit_simpson)[1],coef(fit_shannon)[1],coef(fit_bp)[1])
    #create a vector of values of the slope
    Slope <- rbind(Slope, coef(fit_richness)[2],coef(fit_chao1)[2],coef(fit_simpson)[2],coef(fit_shannon)[2],coef(fit_bp)[2])
    #create a vector of values of the R2
    R2 <- rbind(R2, summary(fit_richness)$r.squared,summary(fit_chao1)$r.squared,summary(fit_simpson)$r.squared,summary(fit_shannon)$r.squared,summary(fit_bp)$r.squared)
    #create a vector of values of the p-values (need to do anova to get this)
    #https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
    P.val <- rbind(P.val,anova(fit_richness)$'Pr(>F)'[1],anova(fit_chao1)$'Pr(>F)'[1],anova(fit_simpson)$'Pr(>F)'[1],anova(fit_shannon)$'Pr(>F)'[1],anova(fit_bp)$'Pr(>F)'[1])
  }



#create a matrix of the the linear model results
lmresults  <-  as.data.frame(cbind(lmresults,Intercept, Slope, R2, P.val))
names(lmresults) <- c("Intercept","Slope","R^2","P.value")
lmresults$Microbe <- c(rep("Bacteria",15),rep("Fungi",15))
lmresults$Timepoint <-c(rep(c("T1"),5),rep(c("T2"),5),rep(c("T3"),5))
lmresults$RichnessIndex <- rep(c("Richness","Chao1","Simpson","Shannon","Berger-Parker"),6)
#now look at my results with names
lmresults

#do an example to check
fit_bac_T3_shannon  <- lm(bacrichness_T3$shannon_e ~ bacrichness_T3$zOTU_shannon_e)
summary(fit_bac_T3_shannon )
fit_fungi_T3_richness  <- lm(fungirichness_T3$richness ~ fungirichness_T3$zOTU_richness)
summary(fit_fungi_T3_richness)

#write csv to a file
write.csv(lmresults, "results/alphadiversitycomparisons_otuvzotu_fungiandbacteria.csv")

########################################################################################################################
######Make loop to get Pearson correlations
#########################################################################################################################
#create empty vectors to store the data
cor_results <- NULL
Pearson.r <- NULL
Pearson.pval<- NULL


#run for loop to do Pearson correlations
for(k in 1:length(mylist_otu)){
  #do the linear model test of otu richness against zotu richness
  cor_richness <- cor.test(mylist_otu[[k]][[9]],mylist_otu[[k]][[24]])
  #do the linear model test of otu richness against zotu chao1
  cor_chao1<- cor.test(mylist_otu[[k]][[3]],mylist_otu[[k]][[18]])
  #do the linear model test of otu richness against zotu Simpson
  cor_simpson <- cor.test(mylist_otu[[k]][[11]],mylist_otu[[k]][[26]])
  #do the linear model test of otu richness against zotu Shannon
  cor_shannon <- cor.test(mylist_otu[[k]][[12]],mylist_otu[[k]][[27]])
  #do the linear model test of otu richness against zotu berger parker
  cor_bp <- cor.test(mylist_otu[[k]][[1]],mylist_otu[[k]][[16]])
  #create a vector of values of the Pearson correlation
  Pearson.r <- rbind(Pearson.r, cor_richness$estimate,cor_chao1$estimate, cor_simpson$estimate, cor_shannon$estimate,cor_bp$estimate )
  #create a vector of values of the pvalues
  Pearson.pval<- rbind(Pearson.pval,cor_richness$p.value,cor_chao1$p.value, cor_simpson$p.value, cor_shannon$p.value,cor_bp$p.value)
}


#create a matrix of the the linear model results
cor_results  <-  as.data.frame(cbind(cor_results,Pearson.r,Pearson.pval))
names(cor_results) <- c("Pearson r","P.value")
cor_results$Microbe <- c(rep("Bacteria",15),rep("Fungi",15))
cor_results$Timepoint <-c(rep(c("T1"),5),rep(c("T2"),5),rep(c("T3"),5))
cor_results$RichnessIndex <- rep(c("Richness","Chao1","Simpson","Shannon","Berger-Parker"),6)
#now look at my results with names
cor_results

#do an example to check
#do an example to check
cortest1 <- cor.test(bacrichness_T3$shannon_e,bacrichness_T3$zOTU_shannon_e)
cortest1$estimate
cortest1$p.value

cortest2 <- cor.test(fungirichness_T3$richness,fungirichness_T3$zOTU_richness)
cortest2$estimate
cortest2$p.value

cor.test(fungirichness_T1$chao1,fungirichness_T1$zOTU_chao1)

#write csv to a file
write.csv(cor_results,"results/Pearsoncorrelations_otuvzotu_fungiandbacteria.csv")

#can repeat with it not rarefied, with filtered OTU tables, with only 20bp removed
###############################################################################################################
######Make total distribution of Pearson r correlations figures
#####################################################################################################################

#correlations <- read.csv("data/correlationsforotuvzotu.csv", header=TRUE)


correlationhistogram <- ggplot(data=cor_results, aes(Pearson.r))+
  geom_histogram(bins=100) +
  scale_x_continuous(limits = c(0.85,1 ))+
  theme(strip.text.x = element_text(size = 14, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text.x=element_text(size=12, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=12),  #make y axis tick sizes bigger
        axis.title=element_text(size=14)) #make y axis label larger

correlationhistogram + facet_wrap(~Microbe)

correlationboxplot <- ggplot(data=cor_results, aes(x=Microbe, y=Pearson.r, col=Microbe))+
  geom_boxplot()+  scale_y_continuous(limits = c(0.5,1 ))+
  labs(x=" ",y="Pearson r", title="C") + # change labels
  theme(axis.text.x=element_text(size=14, angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        axis.title=element_text(size=14),
        plot.title = element_text( face="bold", size=20, hjust=0)) #make y axis label large

correlationboxplot 

correlationboxplot_bw <- ggplot(data=cor_results, aes(x=Microbe, y=Pearson.r))+
  geom_boxplot()+  scale_y_continuous(limits = c(0.5,1 ))+
  labs(x="Microbe",y="Pearson r", title="C") + # change labels
  theme_bw()+
  theme(axis.text.x=element_text(size=14, angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        axis.title=element_text(size=14),
        plot.title = element_text( face="bold", size=20, hjust=0)) #make y axis label large


correlationboxplot_bw 
#base R
boxplot(correlations$Pearson.r ~ correlations$Microbe, ylim=c(0.5,1))

###############################################################################################################
######Combine fungi and bacteria and subset just T1 richness 
#####################################################################################################################
#make facet plot for richness
#add bacteria and fungi names
bacrichness$Microbe <- rep("Bacteria", nrow(bacrichness))
fungalrichness$Microbe <- rep("Fungi", nrow(fungalrichness))

#combine
richness_all <- rbind(bacrichness, fungalrichness)

#subset relevant time points and values
names(richness_all)
myvars <- c("richness", "Timepoint", "Microbe", "zOTU_richness")
richness_all_2 <- richness_all[myvars]
richness_all_T1 <- richness_all_2[richness_all_2$Timepoint=="T1", ]

#make figure
richness <- ggplot(richness_all_T1, aes(x=richness, y=zOTU_richness)) + geom_point(shape=1) +
  geom_smooth(method = "lm") + #add linear smoother
  labs(x="OTU richness",y="ESV richness") + #change x and y label axes
  theme(strip.text.x = element_text(size = 12)) +#make facet labels bigger 
  geom_abline(intercept = 0, slope = 1) #add in 1:1 line


# Divide by time points, going horizontally and wrapping with 3 columns
r <- richness + facet_wrap( ~ Microbe, ncol=2) 

r

###############################################################################################################
######Make just bacteria figure
#####################################################################################################################
###make a bacterial richness figure = panel A
bacrichness_T1 <- bacrichness[bacrichness$Timepoint=="T1", ]

richness_bac_T1 <- ggplot(bacrichness_T1, aes(x=richness, y=zOTU_richness)) + geom_point(shape=1) +
  geom_smooth(method = "lm") + #add linear smoother
  theme_bw()+
  stat_cor(method="pearson", label.x=200, label.y.npc=0.25)  +
  labs(x="Bacterial OTU",y="Bacterial ESV",title="A") + #change x and y label axes
  theme(axis.title=element_text(size=14), axis.text.x=element_text(size=14, angle=70, hjust=1), axis.text.y=element_text(size=14), #make labels bigger 
  plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger 
  
richness_bac_T1


richness_bac <- ggplot(bacrichness, aes(x=richness, y=zOTU_richness, group=Timepoint, color=Timepoint)) + geom_point(shape=1) +
  geom_smooth(method = "lm") + #add linear smoother
  theme_bw()+
  labs(x="Bacterial OTU",y="Bacterial ESV",title="A") + #change x and y label axes
  theme(axis.title=element_text(size=14),axis.text.x=element_text(size=14, angle=70, hjust=1), axis.text.y=element_text(size=14),
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger 


richness_bac

levels(bacrichness$Timepoint) <- c("6 months","12 months","18 months")
bacscatter <- ggscatter(bacrichness, x = "richness", y = "zOTU_richness",
                           add = "reg.line",                         # Add regression line
                           conf.int = TRUE,                          # Add confidence interval
                           color = "Timepoint",            # Color by groups "cyl"
                           shape = "Timepoint"                             # Change point shape by groups "cyl"
)+
  #stat_cor(aes(color = Timepoint), label.x = 200, label.y.npc=0.1, method="pearson")     + # Add correlation coefficient
  labs(x="Bacterial OTU Richnesss",y="Bacterial ESV Richness", title="A") + # change labels
  theme(axis.text.x=element_text(size=14, angle=70, hjust=1), axis.text.y=element_text(size=14), #make y axis tick sizes bigger
        axis.title=element_text(size=14),
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger ) #make y axis label larger

bacscatter
###############################################################################################################
######Make just fungifigure
#####################################################################################################################
###make a bacterial richness figure = panel A
fungalrichness_T1 <- fungalrichness[fungalrichness$Timepoint=="T1", ]

richness_fungi_T1 <- ggplot(fungalrichness_T1, aes(x=richness, y=zOTU_richness)) + geom_point(shape=1) +
  geom_smooth(method = "lm") + #add linear smoother
  theme_bw()+
  stat_cor(method="pearson", label.x=160, label.y.npc=0.25)  +
  labs(x="Fungal OTU",y="Fungal ESV", title="B") + #change x and y label axes
  theme(axis.title=element_text(size=14),axis.text.x=element_text(size=14, angle=70, hjust=1), axis.text.y=element_text(size=14),
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger 
richness_fungi_T1


richness_fungi <- ggplot(fungalrichness, aes(x=richness, y=zOTU_richness, group=Timepoint, color=Timepoint)) + geom_point(shape=1) +
  geom_smooth(method = "lm") + #add linear smoother
  theme_bw()+
  labs(x="Fungal OTU",y="Fungal ESV", title="B") + #change x and y label axes
  theme(axis.title=element_text(size=14), axis.text.x=element_text(size=14, angle=70, hjust=1), axis.text.y=element_text(size=14),
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger  

richness_fungi
#figure out what colors ggplot2 used to build this figure
ggplot_build(richness_fungi)$data[[1]][[1]]

#useful webpage: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/#comments-list

#change the levels from T1 T2 T3 to 6,12,18 months
levels(fungalrichness$Timepoint) <- c("6 months","12 months","18 months")

fungalscatter <- ggscatter(fungalrichness, x = "richness", y = "zOTU_richness",
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          color = "Timepoint",            # Color by groups "cyl"
          shape = "Timepoint"                             # Change point shape by groups "cyl"
)+
  #stat_cor(aes(color = Timepoint),label.x = 100,label.y.npc=0.1, method="pearson")   + # Add correlation coefficient
  labs(x="Fungal OTU Richness",y="Fungal ESV Richness", title="B") + # change labels
  theme(axis.title=element_text(size=14), axis.text.x=element_text(size=14, angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        plot.title = element_text( face="bold", size=20, hjust=0)) #make labels bigger  

fungalscatter







###############################################################################################################
######Make panel
#####################################################################################################################
#####put them all in the same pdf
#really helpful: http://www.sthda.com/english/rpkgs/ggpubr/index.html
#install.packages("ggpubr")
library(ggpubr)


figure1 <- ggarrange(richness_bac, richness_fungi,correlationboxplot_bw , ncol=3, common.legend=TRUE,legend="right")
figure2 <- ggarrange(bacscatter, fungalscatter,correlationboxplot_bw, ncol = 3, common.legend=TRUE,legend="right")
figure3 <- ggarrange(richness_bac_T1,richness_fungi_T1,correlationboxplot_bw, ncol = 3, common.legend=TRUE,legend="right")

figure4 <- ggarrange(bacscatter, fungalscatter, ncol = 2, common.legend=TRUE,legend="bottom")


figure1
figure2
figure3
figure4 

figure1 %>%
  ggexport(filename = "Figures/zotuvotu_richness/Figure1_alpharichness_panel_option1.pdf")


figure2 %>%
  ggexport(filename = "Figures/zotuvotu_richness/Figure1_alpharichness_panel_option2.pdf")

figure3 %>%
  ggexport(filename = "Figures/zotuvotu_richness/Figure1_alpharichness_panel_option3.pdf")

pdf("Figures/zotuvotu_richness/Figure1_alpharichness_panel_option4.pdf", height=6, width=10)
figure4 
dev.off()
