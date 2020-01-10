#!/bin/R
#
# AUTHOR: Olena Maiakovska
#
# Date: 1.09.2019
#
# Object: Selection of polymorphic sites, generation of re-coded matrix of variants, building phylognetic tree and PCA analysis 
#
#LIBRARIES:
library(ape)
library(phytools)
library(ggplot2)
library(ggfortify)
library(pheatmap)

# Dataset  
D_Table <-read.table("~/SNP.vcf", header = F)

# Re-coding the genotypes to alternative allele dosage:
Re_Coded <- apply(D_Table, 2, function(i) {
  i <- sub("^(0/0/0).*", "0", i)
  i <- sub("^(0/0/1).*", "1", i)
  i <- sub("^(0/1/1).*", "2", i)
  i<- sub("^(1/1/1).*", "3", i)})

#Phylogenetic tree
Re_coded = Hetero
snp = rbind(HD1=as.numeric(Hetero$HD1), Moosweiher=as.numeric(Hetero$Moosweiher), madAla=as.numeric(Hetero$madAla), madAna=as.numeric(Hetero$madAna), madIta=as.numeric(Hetero$madIta), madVak=as.numeric(Hetero$madVak), Estonia = as.numeric(Hetero$Estonia), Reilingen=as.numeric(Hetero$Reilingen),Czech=as.numeric(Hetero$Czech_Vodnany), Hungary=as.numeric(Hetero$Hungary), Malta=as.numeric(Hetero$Malta), Romania=as.numeric(Hetero$Romania), Ukraine=as.numeric(Hetero$Ukraine), Singliser_See=as.numeric(Hetero$Singliser_See), Slovakia=as.numeric(Hetero$Slovakia))
stree <- nj(dist.gene(snp))
plot(unroot(stree),type="unrooted",no.margin=TRUE, edge.width=2, show.tip.label = T, use.edge.length = F)
add.scale.bar(x=5, y=-300 , cex = 0.7, font = 1, col="red")


#PCA analysis and Heatmap clusterization:
Names <- read.csv("~/PopulationsNames.csv", sep= " ", header=F)
rownames(Names) <- Names
df <- t(Re_Coded)  ### data: rows - Names and columns - SNPs
tdata <- cbind(data, df)
autoplot(prcomp(df), data=tdata, label.size = 4,   label = TRUE) 

#Heatmap:
metadata_populations <- data.frame(data$V2)
rownames(metadata_populations) = data$V1
pheatmap(DataFile, cluster_rows = T, cluster_cols = T, kmeans_k = 600, iter.max = 100, annotation_col = metadata_populations,  color = c("lightblue4","lightblue1","darkred"))






