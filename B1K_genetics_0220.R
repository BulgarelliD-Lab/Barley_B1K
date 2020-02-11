#############################################################
#
# Ref to the ARTICLE
# 
# Rodrigo Alegria Terrazas (r.z.alegriaterrazas@dundee.ac.uk) and Davide Bulgarelli (d.bulgarelli@dundee.ac.uk)
# 11/02/20
# script to reproduce calculations and figures presented in the manuscript
# 
# 
#############################################################

#############################################################
# Libraries and functions required
#############################################################

#These initial commands are required to clean-up the memory and start a new session
rm(list=ls())
dev.off()

#set working directory
#setwd("C:/Users/Ra42320/Box/Davide_lab_manuscripts/Rodrigo_B1K_2019/R_data")
setwd("C:/Users/DB42008/Box Sync/Davide_lab_manuscripts/Rodrigo_B1K_2019/R_data/")
#laptop
#setwd("C:/Users/Ra42320/Box Sync/Davide_lab_manuscripts/Rodrigo_B1K_2019/R_data")

#library location
.libPaths()

#load the required packages 
library("phyloseq")
library("ggplot2")
library("vegan")
library("ade4")
library("adegenet")
library("harrietr")
library("ggpubr")

#R session info
sessionInfo()

#import the Phyloseq file (see R code phyloseq)
JH04_data_phyloseq_prop <- readRDS(file = "JH04_data_phyloseq_prop.rds")
JH04_data_phyloseq_prop

#extract OTU table and 
otu_table_prop <- as.data.frame(otu_table(JH04_data_phyloseq_prop))
otu_table_prop[1:10, ]

#extract the design file
design_prop <- as.data.frame(sample_data(JH04_data_phyloseq_prop))

#inspect the file
colnames(otu_table_prop)
rownames(design_prop)

#replace column names in the OTU table
colnames(otu_table_prop) <- design_prop$B1K
otu_table_prop[1:10, ]

#create a new otu table that combines the average values of the individual samples
otu_table_prop_average <- as.data.frame(do.call(cbind, by(t(otu_table_prop),INDICES=names(otu_table_prop),FUN=colMeans)))

#import a new design file where row names are the individual B1K identifiers
design_average <- read.delim("Map_JH04_24_B1K_ids.txt", sep = "\t", header=TRUE, row.names=1)
JH04_average_map <- sample_data(design_average)

#subset the taxonomy file for the OTU present in the dataset
average_taxa_ordered <- tax_table(JH04_data_phyloseq_prop)

#extract the rooted phylogentic tree
tree_average = phy_tree(JH04_data_phyloseq_prop)

#create a new phyloseq object
#a) The OTU Table counts
JH04_OTU_average <- otu_table(otu_table_prop_average, taxa_are_rows=TRUE)

#merge the files and create the phyloseq object
JH04_average_data_phyloseq <- merge_phyloseq(JH04_OTU_average, average_taxa_ordered, JH04_average_map,  tree_average)

#import the genetic distnace
#this file is named Supplementary worksheet ws 10 in Additional file 2
GD_data <- read.delim("simple_matching_GD_all_2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
class(GD_data)


#replace the NA values
GD_data[is.na(GD_data)] <- 0
#change the sample ID
#this file is named Supplementary worksheet ws 11 in Additional file 2
design_average_GD_24 <- read.delim("Map_JH04_24_B1K_ids.txt",sep = "\t", header=TRUE, row.names=1 )

#intersect with 24 samples
Microbe_GD_24 <- intersect(rownames(design_average_GD_24), rownames(GD_data))
#write(Microbe_GD_24, file="Microbe_GD_24.txt")
colnames_GD_Data<- colnames(GD_data)
#write(colnames_GD_Data, file="colnames_GD_Data.txt")
Microbe_GD_24_2 <- c("Barke", "Morex", "Steptoe", "Bowman", "B1K.02.18", "B1K.03.09", "B1K.04.04", "B1K.05.13", "B1K.08.18", "B1K.11.11", "B1K.12.10",
                     "B1K.14.04","B1K.15.19", "B1K.17.10", "B1K.18.16", "B1K.20.13", "B1K.21.11", "B1K.30.07", "B1K.31.01", "B1K.33.03", "B1K.34.20", "B1K.35.11",
                     "B1K.37.06", "B1K.48.06")
GD_data_24 <- GD_data[Microbe_GD_24_2, Microbe_GD_24_2]

#create a distance matrix
GD_24_dist <- as.dist(GD_data_24)

#cluster dendrogram (Additional file 1: Figure S4)
GD_24_clust <- hclust(GD_24_dist, method = "average")
plot(GD_24_clust)

#BC distance
BC <- phyloseq::distance(JH04_average_data_phyloseq, "bray")
WU <- phyloseq::distance(JH04_average_data_phyloseq, "unifrac", weighted= TRUE)

#Subset wild
Wild_samples <- rownames(design_average_GD_24)[which(design_average_GD_24$Description != "Modern")]
design_average_GD_20 <- design_average_GD_24[Wild_samples, ]
#intersect with 20 samples
Microbe_GD_20 <- intersect(rownames(design_average_GD_20), rownames(GD_data))
GD_data_20 <- GD_data[Microbe_GD_20, Microbe_GD_20]
GD_20_dist <- as.dist(GD_data_20)

#subset the BC distance for wild barleys
JH04_average_data_phyloseq_rhizo_20 <- prune_samples(Microbe_GD_20, JH04_average_data_phyloseq)
BC_20 <- phyloseq::distance(JH04_average_data_phyloseq_rhizo_20, "bray")
WU_20 <- phyloseq::distance(JH04_average_data_phyloseq_rhizo_20, "unifrac", weighted= TRUE)

#subset modern
M_samples <- rownames(design_average_GD_24)[which(design_average_GD_24$Geo== "M")]
design_M <- design_average_GD_24[M_samples, ]
Microbe_GD_M <- intersect(rownames(design_M), rownames(GD_data_24))
GD_data_M <- GD_data[Microbe_GD_M,Microbe_GD_M ]
GD_M_dist <- as.dist(GD_data_M)

#subset M
JH04_average_data_phyloseq_rhizo_M <- subset_samples(JH04_average_data_phyloseq, Geo == "M")
BC_M <- phyloseq::distance(JH04_average_data_phyloseq_rhizo_M, "bray")
WU_M <- phyloseq::distance(JH04_average_data_phyloseq_rhizo_M,"unifrac", weighted= TRUE)


#####mantel test

#All samples
#BC distance
mantel.rtest(BC, GD_24_dist, nrepet = 9999)
#WU distance
mantel.rtest(WU, GD_24_dist, nrepet = 9999)

#Modern
mantel.rtest(GD_M_dist, BC_M, nrepet = 9999)
#WU distance
mantel.rtest(GD_M_dist, WU_M, nrepet = 9999)

#Only wild
mantel.rtest(GD_20_dist, BC_20 , nrepet = 9999)
#wu distance
mantel.rtest(GD_20_dist, WU_20, nrepet = 9999)

#plotting
WU_20_vector <- melt_dist(as.matrix(WU_20), order = NULL, dist_name = "WU_dist")
B1K_GD_20_vector <- melt_dist(as.matrix(GD_20_dist), order = NULL, dist_name = "GD_dist")
WU_20_vector
B1K_GD_20_vector

#merge the datasets
Correlation_data_WU <- cbind(WU_20_vector[, "WU_dist"], B1K_GD_20_vector[, "GD_dist"])
colnames(Correlation_data_WU) <- c("WU_dist", "GD_dist")
Correlation_data_WU <- as.data.frame(Correlation_data_WU)

#plotting
dev.off()
ggscatter(Correlation_data_WU, x ="GD_dist" , y ="WU_dist" , add = "reg.line", conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue", fill = "lightgray"))

#end of the code




