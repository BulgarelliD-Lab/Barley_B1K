#############################################################
#
# Ref to the ARTICLE
# 
# Rodrigo Alegria Terrazas (r.z.alegriaterrazas@dundee.ac.uk) and Davide Bulgarelli (d.bulgarelli@dundee.ac.uk)
# 22/04/20
# script to reproduce calculations and figures presented in the manuscript
# 
# 
###########################################################
#ANCOM 
#Analysis for revisions of SciRep 04/20
###########################################################

rm(list=ls())
dev.off()

#set the working directory
setwd("/cluster/db/ralegriaterrazas/R_data_JH04")

#############################################################
# Libraries required
#############################################################

library(exactRankTests)
library(nlme)

#ANCOM package
#https://github.com/sidhujyatha/ANCOM

source("ANCOM_updated_code.R", echo=F)

#retrieve R and package versions 
sessionInfo()



# Step 1:  Input files

#This table was generated in QIIME w/o plant-derived OTUs and singletons with absolute number of reads at phylum level
#file "Map_JH04_24_L2.txt"  is  worksheet 4 in Supplementary Dataset 1
#file "Map_ANCOM.txt" reduced version of mapping file for main analysis 
#Samples names need to be written as "Sample.ID" is worksheet 5 in Supplementary Dataset 1

Phylum_table<- read.delim("Map_JH04_24_L2.txt", sep = "\t", header=TRUE)

design_ANCOM <-read.delim("Map_ANCOM.txt", sep = "\t", header=TRUE)



# Step 2: ANCOM

comparison_test=ANCOM.main (OTUdat=Phylum_table,
                            Vardat=design_ANCOM[design_ANCOM$Microhabitat%in%c("Rhizosphere","Bulk"), ],
                            adjusted=F, # Not adjusting for anything
                            repeated=F, # No repeated data
                            main.var="Microhabitat", #Variable Microhabitat
                            adj.formula=NULL, # Not adjusting for anything
                            repeat.var=NULL, # No repeated data
                            longitudinal=FALSE, # Just 1 timepoint
                            random.formula=NULL,
                            multcorr=2, # Taxa-based p-value correction (recommended)
                            sig=0.05, # $\alpha$ = 0.05
                            prev.cut=0.90) # Remove taxa absent from >=90% of samples

W.taxa_Microhabitat = comparison_test$W.taxa
write.table(W.taxa_Microhabitat, "ANCOM2_Microhabitat.txt", sep="\t")
# the above output is saved as worksheet 6 in Supplementary Dataset 1





