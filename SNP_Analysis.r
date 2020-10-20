#!/bin/R
#
# AUTHOR: Olena Maiakovska
#
# Date: 1.09.2019
#
# Object: Filtering of heterozygous positions, selection of polymorphic sites, generation of re-coded matrix of variants, building phylognetic tree and PCA analysis 
#
#LIBRARIES:
require(reshape)
library(ape)
library(phytools)
library(ggplot2)
library(ggfortify)
library(pheatmap)

#vcf files consisting either European populations +reference HD1 or European + Malagasy tohether or all populations together with Reilingen samples
# Dataset vcf file with P. virginalis samples including reference HD1  
VCF_data <-read.table("~/Multisample.vcf", header = F)
#Selection of genotypes from the vcf table:
VCF_Data <- VCF_Data[,10:22]
#Filtering for heterozygous positions in reference genotype:
Hetero <- VCF_Data[grep("0/0/0", VCF_Data$HD1),]

# Re-coding the genotypes to alternative allele dosage:
Re_Coded <- apply(Hetero, 2, function(i) {
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


#LOH analysis:

Czech_Vodnany = colsplit(Genotypes[,1], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Estonia = colsplit(Genotypes[,2], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
HD1 = colsplit(Genotypes[,3], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Hungary = colsplit(Genotypes[,4], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Malta = colsplit(Genotypes[,5], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Moosweiher = colsplit(Genotypes[,6], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Reilingen = colsplit(Genotypes[,7], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Romania = colsplit(Genotypes[,8], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Singliser_See = colsplit(Genotypes[,9], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Slovakia = colsplit(Genotypes[,10], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))
Ukraine=colsplit(Genotypes[,11], split = "\\:", names = c('GT', 'AO','DP', 'QA','QR','RO'))

Common=cbind(Czech_Vodnany,Estonia, HD1, Hungary,Malta,Moosweiher,Reilingen,Romania,Singliser_See,Slovakia,Ukraine)
names(Common) <- c('GT', 'AO','DP', 'QA','QR','RO', 
                   'GT', 'AO1','DP1','QA','QR','RO1',
                   'GT', 'AO2','DP2', 'QA','QR','RO2',
                   'GT', 'AO3','DP3', 'QA','QR','RO3',
                   'GT', 'AO4','DP4', 'QA','QR','RO4',
                   'GT', 'AO5','DP5', 'QA','QR','RO5',
                   'GT', 'AO6','DP6', 'QA','QR','RO6',
                   'GT', 'AO7','DP7', 'QA','QR','RO7',
                   'GT', 'AO8','DP8', 'QA','QR','RO8',
                   'GT', 'AO9','DP9', 'QA','QR','RO9',
                   'GT', 'AO10','DP10','QA','QR','RO10')



Common =subset(Common, RO>5)
Common=subset(Common, RO1>5)
Common=subset(Common, RO2>5)
Common=subset(Common, RO3>5)
Common=subset(Common, RO4>5)
Common=subset(Common, RO5>5)
Common=subset(Common, RO6>5)
Common=subset(Common, RO7>5)
Common=subset(Common, RO8>5)
Common=subset(Common, RO9>5)
Common=subset(Common, RO10>5)

VAF_Czech_Vodnany = Common$AO/Common$DP
VAF_Estonia = Common$AO1/Common$DP1
VAF_HD1 = Common$AO2/Common$DP2
VAF_Hungary =  Common$AO3/Common$DP3
VAF_Malta = Common$AO4/Common$DP4
VAF_Moosweiher= Common$AO5/Common$DP5
VAF_Reilingen = Common$AO6/Common$DP6
VAF_Romania = Common$AO7/Common$DP7
VAF_Singliser_See = Common$AO8/Common$DP8
VAF_Slovakia = Common$AO9/Common$DP9
VAF_Ukraine = Common$AO10/Common$DP10


seq <- seq(1, length(LOHCz), 1000 )

LOHCz <-  (VAF_Czech_Vodnany)
LOHCz_Vector <- sapply(seq, function(i) {mean(LOHCz[i:(i+1000)])})


LOHEst <-VAF_Estonia
LOHEst_Vector <- sapply(seq, function(i) {mean(LOHEst[i:(i+1000)])})


LOHHD1 <-  (VAF_HD1)
LOHHD1_Vector <- sapply(seq, function(i) {mean(LOHHD1[i:(i+1000)])})

LOHMoosweiher <- VAF_Moosweiher
LOHMoosweiher_Vector <- sapply(seq, function(i) {mean(LOHMoosweiher[i:(i+1000)])})

LOHMalta <- VAF_Malta
LOHMalta_Vector <- sapply(seq, function(i) {mean(LOHMalta[i:(i+1000)])})

LOHHungary <- VAF_Hungary
LOHHungary_Vector <- sapply(seq, function(i) {mean(LOHHungary[i:(i+1000)])})


LOHReilingen <- VAF_Reilingen
LOHReilingen_Vector <- sapply(seq, function(i) {mean(LOHReilingen[i:(i+1000)])})

LOH_Romania <-  VAF_Romania
LOHRomania_Vector <- sapply(seq, function(i) {mean(LOH_Romania[i:(i+1000)])})

LOH_Singliser_See <- VAF_Singliser_See
LOHSingliser_See_Vector <- sapply(seq, function(i) {mean(LOH_Singliser_See[i:(i+1000)])})

LOH_Slovakia <- VAF_Slovakia
LOH_Slovakia_Vector <- sapply(seq, function(i) {mean(LOH_Slovakia[i:(i+1000)])})

LOH_Ukraine <-  VAF_Ukraine
LOH_Ukraine_Vector <- sapply(seq, function(i) {mean(LOH_Ukraine[i:(i+1000)])})

plot(LOHHD1_Vector, type="l", lwd = 1,col = "red", ylab = "VAF", xlab = "kb", ylim = c(0.1,0.5))
lines(LOHSingliser_See_Vector, col = "black",lwd = 0.3,)
lines(LOHMoosweiher_Vector, col="black",lwd = 0.3,)
lines(LOHMalta_Vector, col = "black",lwd = 0.3,)
lines(LOHHungary_Vector, col = "black",lwd = 0.3,)
lines(LOHReilingen_Vector, col = "black",lwd = 0.3,)
lines(LOHRomania_Vector, col = "black",lwd = 0.3,)
lines(LOHSingliser_See_Vector, col = "black",lwd = 0.3,)
lines(LOH_Slovakia_Vector, col = "black",lwd = 0.3,)
lines(LOH_Ukraine_Vector, col = "black",lwd = 0.3,)

