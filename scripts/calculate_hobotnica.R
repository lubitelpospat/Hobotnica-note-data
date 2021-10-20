#!/usr/bin/env Rscript

library(Hobotnica)
library(edgeR)
args = commandArgs(trailingOnly=TRUE)


distMatrix <- args[2]

annotationFile <- args[1]



anno_raw = read.table(annotationFile, header=TRUE, sep=",")
#print("Processing...")
colnames(anno_raw) <- c("Run", "source_name")
anno = data.frame(anno_raw$source_name, row.names = anno_raw$Run)
colnames(anno) <- "group"

cm = read.table(distMatrix, header=TRUE, sep=",")
cm = as.matrix(data.frame(cm[, -1], row.names = cm[, 1]))

colnames(cm) <- rownames(cm)
#cm<- cm[rownames(anno), rownames(anno)]
print(dim(cm))
#cm_cpm <- cpm(cm)
print(Hobotnica(cm, anno$group))


