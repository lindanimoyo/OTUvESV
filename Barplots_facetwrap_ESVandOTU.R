#March 6, 2018
#Make facet wrap figures for bacterial otu v esv top genus

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

#################################################################
#############read in 16S esv vs OTU summarized by GENUS   
#################################################################
#read in bacteria 16S otu versus esv by genera        
genus_sorted_16S_otu <- read.csv("Figures/otuvzotu_taxonomy/16S_Inoculum_phylasorted_genus_otu_rarefied.csv", row.names=1)
genus_sorted_16S_esv <- read.csv("Figures/otuvzotu_taxonomy/16S_Inoculum_phylasorted_genus_zotu_rarefied.csv", row.names=1)

#bind them together to a single data frame
genus_sorted_16S_all <- rbind(genus_sorted_16S_otu,genus_sorted_16S_esv) 

#re order sites
genus_sorted_16S_all$Sitefullnames <- factor(genus_sorted_16S_all$Sitefullnames, 
                                           levels = c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))  

################################################
#############plot figure
##############################################
genus_barplot1_otuvzotu<-
  ggplot(genus_sorted_16S_all, aes(x=Sitefullnames, y=mean, fill=Phylum)) + 
  geom_bar(stat="identity",color="black", width=0.4) +
  scale_fill_brewer(palette = "Paired")+ 
  #When you use scale_fill_brewer(palette = "Paired"), the maximum is 12, or 13 ("other phylum" is not colored)
  labs(x="Site", y="Bacteria Relative abundance (%)", fill="Genus") + theme_bw()+
  theme(strip.text = element_text(size=14),#make facet labels larger
    axis.text.x=element_text(size=12,angle=30, hjust=1),
        axis.text.y=element_text(size=14),
        legend.text = element_text(size=12)
  )

plot_genus_bacteria <- genus_barplot1_otuvzotu + facet_wrap(~otuvesv)

plot_genus_bacteria

ggsave(plot_genus_bacteria,filename=paste("Figures/otuvzotu_taxonomy/genusbarplot_inoculum_bacteria_otuvzotu.pdf"), width=10, height=5)


#################################################################
#############read in 16S esv vs OTU summarized by FAMILY
#################################################################
#read in bacteria 16S otu versus esv by family       
fam_sorted_16S_otu <- read.csv("Figures/otuvzotu_taxonomy/16S_Inoculum_phylasorted_rarefied_otu.csv", row.names=1)
fam_sorted_16S_esv <- read.csv("Figures/otuvzotu_taxonomy/16S_Inoculum_phylasorted_zotu_rarefied.csv", row.names=1)

#bind them together to a single data frame
fam_sorted_16S_all <- rbind(fam_sorted_16S_otu,fam_sorted_16S_esv) 

#re order sites
fam_sorted_16S_all$Sitefullnames <- factor(fam_sorted_16S_all$Sitefullnames, 
                                               levels = c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))   
################################################
#############plot figure
#############################################

fam_barplot1_otuvzotu<-
  ggplot(fam_sorted_16S_all , aes(x=Sitefullnames, y=mean, fill=Phylum)) + 
  geom_bar(stat="identity",color="black", width=0.4) +
  scale_fill_brewer(palette = "Paired")+ 
  #When you use scale_fill_brewer(palette = "Paired"), the maximum is 12, or 13 ("other phylum" is not colored)
  labs(x="Site", y="Bacteria Relative abundance (%)", fill="Family") + theme_bw()+
  theme(strip.text = element_text(size=14),#make facet labels larger
        axis.text.x=element_text(size=12,angle=30, hjust=1),
        axis.text.y=element_text(size=14),
        legend.text = element_text(size=12)
  )

plot_family_bacteria <-fam_barplot1_otuvzotu + facet_wrap(~otuvesv)
plot_family_bacteria 

ggsave(plot_family_bacteria ,filename=paste("Figures/otuvzotu_taxonomy/fambarplot_inoculum_bacteria_otuvzotu.pdf"), width=10, height=5)


#################################################################
#############read in ITS2 fungi esv vs OTU summarized by GENUS
#################################################################
#read in bacteria 16S otu versus esv by genera        
genus_sorted_fungi_otu <- read.csv("Figures/otuvzotu_taxonomy/Fungi_Inoculum_genus_phylasorted_otu.csv", row.names=1)
genus_sorted_fungi_esv <- read.csv("Figures/otuvzotu_taxonomy/Fungi_Inoculum_genus_phylasorted_zotu.csv", row.names=1)

#bind them together to a single data frame
genus_sorted_fungi_all <- rbind(genus_sorted_fungi_otu,genus_sorted_fungi_esv) 

#re order sites
genus_sorted_fungi_all$Sitefullnames <- factor(genus_sorted_fungi_all$Sitefullnames, 
                                               levels = c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))   

#sort 
genus_sorted_fungi_all$Phylum <- factor(genus_sorted_fungi_all$Phylum, levels=c( "Dothideaceae;g__Endoconidioma","Dothideomycetes; g__Ramimonilia", "Dothioraceae;g__Aureobasidium" ,"Dothioraceae;g__unidentified",
                                                                                 "Didymellaceae;g__Nothophoma", "Mycosphaerellaceae;g__Mycosphaerella", "Phaeosphaeriaceae;g__Stagonospora", "Pleosporaceae;g__Alternaria",        
                                                                                 "Pleosporales;g__Ascochyta", "Pleosporales;g__Phoma", "Rhytismataceae;g__unidentified",  "Rutstroemiaceae;g__unidentified","other taxa","unidentified Ascomycota"))

################################################
#############plot figure
##############################################
#display.brewer.pal(12, "Paired")
#brewer.pal(12, "Paired")

genus_barplot1_otuvzotu_fungi <-
  ggplot(genus_sorted_fungi_all, aes(x=Sitefullnames, y=mean, fill=Phylum)) + 
  geom_bar(stat="identity",color="black", width=0.4) +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928","darkgrey","white"))+
  #When you use scale_fill_brewer(palette = "Paired"), the maximum is 12, or 13 ("other phylum" is not colored)
  labs(x="Site", y="Fungi Relative abundance (%)", fill="Genus") + theme_bw()+
  theme(strip.text = element_text(size=14),#make facet labels larger
        axis.text.x=element_text(size=12,angle=30, hjust=1),
        axis.text.y=element_text(size=14),
        legend.text = element_text(size=12)
  )

plot_genus_fungi <- genus_barplot1_otuvzotu_fungi + facet_wrap(~otuvesv)

plot_genus_fungi

ggsave(plot_genus_fungi,filename=paste("Figures/otuvzotu_taxonomy/genusbarplot_inoculum_fungi_otuvzotu.pdf"), width=10, height=5)


#################################################################
#############read in ITS2 fungi esv vs OTU summarized by FAMILY
#################################################################
      
fam_sorted_fungi_otu <- read.csv("Figures/otuvzotu_taxonomy/Fungi_Inoculum_familysorted_otu.csv", row.names=1)
fam_sorted_fungi_esv <- read.csv("Figures/otuvzotu_taxonomy/Fungi_Inoculum_familysorted_zotu.csv", row.names=1)


#bind them together to a single data frame
fam_sorted_fungi_all <- rbind(fam_sorted_fungi_otu,fam_sorted_fungi_esv) 

#re order sites
fam_sorted_fungi_all$Sitefullnames <- factor(fam_sorted_fungi_all$Sitefullnames, 
                                               levels = c("Desert","Scrubland","Grassland","Pine-Oak","Subalpine"))     

fam_sorted_fungi_all$Phylum <- factor(fam_sorted_fungi_all$Phylum, levels=c(";f__other families","Ascomycota;f__unidentified",
                                                                           "Ascomycota;f__Dothideaceae","Ascomycota;f__Dothideomycetes","Ascomycota;f__Dothioraceae","Ascomycota;f__Didymellaceae",
                                                                            "Ascomycota;f__Herpotrichiellaceae","Ascomycota;f__Mycosphaerellaceae","Ascomycota;f__Phaeosphaeriaceae",
                                                                            "Ascomycota;f__Pleosporaceae","Ascomycota;f__Pleosporales","Ascomycota;f__Rhytismataceae","Ascomycota;f__Rutstroemiaceae"))
  

################################################
#############plot figure
##############################################
fam_barplot1_otuvzotu_fungi <-
  ggplot(fam_sorted_fungi_all , aes(x=Sitefullnames, y=mean, fill=Phylum)) + 
  geom_bar(stat="identity",color="black", width=0.4) +
  scale_fill_brewer(palette = "Paired")+ 
  #When you use scale_fill_brewer(palette = "Paired"), the maximum is 12, or 13 ("other phylum" is not colored)
  labs(x="Site", y="Fungi Relative abundance (%)", fill="Family") + theme_bw()+
  theme(strip.text = element_text(size=14),#make facet labels larger
        axis.text.x=element_text(size=12,angle=30, hjust=1),
        axis.text.y=element_text(size=14),
        legend.text = element_text(size=12)
  )

plot_fam_fungi <- fam_barplot1_otuvzotu_fungi + facet_wrap(~otuvesv)


plot_fam_fungi

ggsave(plot_fam_fungi,filename=paste("Figures/otuvzotu_taxonomy/fambarplot_inoculum_fungi_otuvzotu.pdf"), width=10, height=5)

