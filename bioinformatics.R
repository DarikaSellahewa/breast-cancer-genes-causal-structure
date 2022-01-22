---
title: "Bioinformatics Project"
author: "Dilinee Sellahewa"
date: "02/10/2021"
---

knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
tinytex::install_tinytex()

library(bnlearn)
library(gRain)
library(pcalg)
library(biclust)
library(dplyr)
library(caTools)
library(e1071)
library(caret)
library(tidyr)
library(klaR)

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

data <- read.csv("BRCA_RNASeqv2_top50.csv",header=T,strip.white=TRUE)
#Removing class variable
gene_data <- subset(data, select = -c(class))

#generate gene network
gene_data <- na.omit(gene_data)
n <- nrow(gene_data) 
V <- colnames(gene_data) #node names
p <- length(V)
suffStat = list(C = cor(gene_data), n= nrow(gene_data))
pc.gene <- pc(suffStat, indepTest = gaussCItest, labels = V,  alpha=0.01)

V.index <-  c()
for (ind in 1: length(V)) {
  V.index <-  c(V.index, ind)
}
#for the representation labels of the graph are numbers and refer label_ref table for respective gene names
pc.di <- pc(suffStat, indepTest = gaussCItest, labels = as.character(V.index),  alpha=0.01)
if (require(Rgraphviz)) {
par(mfrow=c(1,1))
plot(pc.di@graph, main = "Gene Network")
}
label_ref <- cbind(V, V.index)
colnames(label_ref) <- c("Gene name", "Label")
write.csv(label_ref, "label.csv" )

#top 10 genes
BF1_index <- grep("EBF1", colnames(gene_data))
ida_result <- c()
for (v in 1:length(V)) {
  cov_val <- min(abs(ida(v, EBF1_index, cov(gene_data),pc.gene@graph)))
  ida_result <- c(ida_result, cov_val[c(1)])
}
EBF1_genes_index <- tail(order(ida_result),11)
EBF1_top_genes <- c()
for (i in EBF1_genes_index) {
  EBF1_top_genes <- c(EBF1_top_genes,V[i])
}

#markov blanket
M <- learn.mb(gene_data, "ABCA9", method="iamb", alpha=0.01)

#Average of all genes
average_exp <- mean(colMeans(gene_data), na.rm = TRUE)

#cloning the data for convenient use
data_clone <- data

#binarizing data using threshold as the average of all genes.
binary_data <- binarize( gene_data, threshold = average_exp)

data_clone <- data_clone %>% mutate(class = recode(class, "N" = 0, "C" = 1))
 
#pc-simple - parent-childeren set
pc_gene <- pcSelect(data_clone$class, binary_data, alpha=0.05, verbose = 0)
top_genes <- c()
for (i in 1:length(pc_gene$G)) {
  if(pc_gene$G[i] == TRUE)
  { 
    top_genes <- c(top_genes, pc_gene$G[i])
  }
}

#Naive bayes classification using all genes

binary_data <- cbind(binary_data, data_clone$class)#creating the datafrmae
colnames(binary_data)[which(names(binary_data) == "data_clone$class")] <- "classV"

# Define train control for 10 fold cross validation
tc <- trainControl(method="cv", number=10)

#prepare attributes
attr <- subset(binary_data, select = -c(classV))

# Fit Naive Bayes Model
model_1 <- train(attr, as.factor(binary_data$classV), 'nb', trControl=tc)

# Summarise Results
print(model_1)


#Naive bayes classification using parents and children of class variable

pc_genes_index <- tail(order(pc_gene$zMin),length(top_genes))#selecting indices
mat <- matrix(ncol = 1, nrow = nrow(binary_data))#creating empty matrix 
top_binary_data <- data.frame(mat)#creating empty data frame

#filling datafrmae with selected genes
for (g in pc_genes_index) {
  top_binary_data <- cbind(top_binary_data,binary_data[g])
}

top_binary_data <- subset(top_binary_data, select = -c(mat))#removing empty column
top_binary_data <- cbind(top_binary_data, data_clone$class)#adding class variable
colnames(top_binary_data)[which(names(top_binary_data) == "data_clone$class")] <- "classV"
  
# Define train control for 10 fold cross validation
tc <- trainControl(method="cv", number=10)

#prepare attributes
attr <- subset(top_binary_data, select = -c(classV))

# Fit Naive Bayes Model
model_2 <- train(attr, as.factor(top_binary_data$classV), 'nb', trControl=tc)

# Summarise Results
print(model_2)

#bayesian network

cn <- c(1,0)
btnl <- cptable(~btnl, 
                 values=c(sum(binary_data$BTNL9 == 1),sum(binary_data$BTNL9 == 0))
                 ,levels=cn)

cdlg.btnl <- cptable(~cdlg | btnl, 
                 values=c(sum(binary_data$BTNL9 == 1 & binary_data$CD300LG == 1),
                          sum(binary_data$BTNL9 == 1 & binary_data$CD300LG == 0),
                          sum(binary_data$BTNL9 == 0 & binary_data$CD300LG == 1),
                          sum(binary_data$BTNL9 == 0 & binary_data$CD300LG == 0)),
                 levels=cn)


class.cdlg <- cptable(~class | cdlg, 
                 values=c(sum(binary_data$CD300LG == 1 & binary_data$classV == 1),
                          sum(binary_data$CD300LG == 1 & binary_data$classV == 0),
                          sum(binary_data$CD300LG == 0 & binary_data$classV == 1),
                          sum(binary_data$CD300LG == 0 & binary_data$classV == 0)),
                 levels=cn)


igsf.class <- cptable(~igsf | class, 
                 values=c(sum(binary_data$classV == 1 & binary_data$IGSF10 == 1),
                          sum(binary_data$classV == 1 & binary_data$IGSF10 == 0),
                          sum(binary_data$classV == 0 & binary_data$IGSF10 == 1),
                          sum(binary_data$classV == 0 & binary_data$IGSF10 == 0)),
                 levels=cn)

bi_d <- binary_data

abca.igbt <-
  cptable(~abca|igsf:btnl,
          values=c(sum(bi_d$BTNL9 == 1 & bi_d$IGSF10 == 1  & bi_d$ABCA9 == 1),
                   sum(bi_d$BTNL9 == 1 & bi_d$IGSF10 == 1  & bi_d$ABCA9 == 0),
                   sum(bi_d$BTNL9 == 1 & bi_d$IGSF10 == 0  & bi_d$ABCA9 == 1),
                   sum(bi_d$BTNL9 == 1 & bi_d$IGSF10 == 0  & bi_d$ABCA9 == 0),
                   sum(bi_d$BTNL9 == 0 & bi_d$IGSF10 == 1  & bi_d$ABCA9 == 1),
                   sum(bi_d$BTNL9 == 0 & bi_d$IGSF10 == 1  & bi_d$ABCA9 == 0),
                   sum(bi_d$BTNL9 == 0 & bi_d$IGSF10 == 0  & bi_d$ABCA9 == 1),
                   sum(bi_d$BTNL9 == 0 & bi_d$IGSF10 == 0  & bi_d$ABCA9 == 0))
          ,levels=cn)

plist <- compileCPT(list(btnl, cdlg.btnl, class.cdlg, igsf.class, abca.igbt))

#Bayesian network recreated after calculations
net <- grain(plist)
if (require(Rgraphviz)) {
par(mfrow=c(1,1))
plot(net$dag, main = "Bayesian Network")
}

#conditional probabailities
plist$btnl
plist$cdlg
plist$class
plist$igsf
plist$abca

#Probabilities of genes having high expressions
querygrain(net, nodes=c("btnl"), type="marginal")
querygrain(net, nodes=c("cdlg"), type="marginal")
querygrain(net, nodes=c("igsf"), type="marginal")
querygrain(net, nodes=c("abca"), type="marginal")


#To calculate the probability of having a cancer given expression level of CD300LG is high and the expression level of BTNL9 is low
querygrain(net, nodes=c("class", "cdlg", "btnl"), type="conditional")


