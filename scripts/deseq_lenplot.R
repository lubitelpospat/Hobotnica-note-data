library("DESeq2")
library("limma")
library("Hobotnica")

library("doParallel")
library("foreach")
library("amap")

args = commandArgs(trailingOnly=TRUE)

countMatrixFile <- args[1]
annotationFile <- args[2]


cm <- read.table(countMatrixFile, header=T, sep=",")
cm <- data.frame(cm[, -1], row.names=cm[, 1])
print("COUNT MATRIX SHAPE")
print(dim(cm))
anno_raw = read.table(annotationFile, header=TRUE, sep=",")

colnames(anno_raw) <- c("Run", "source_name")
anno = data.frame(anno_raw$source_name, row.names = anno_raw$Run)
colnames(anno) <- "Samples"
print(head(anno))

coldata <- data.frame(condition=as.factor(anno$Samples))
dds <- DESeqDataSetFromMatrix(countData = cm,
                              colData = coldata,
                              design = ~ condition)


cds = estimateSizeFactors(dds);
cds = estimateDispersions(cds);
cds = nbinomWaldTest(cds);

de = data.frame(results(cds));
de_deseq2 = de[!is.na(de$padj), ];


generateScoresByLength <- function(dataset, annotation, all_genes, range=10:200, n_cores=2) {
  
  cl <- parallel::makeForkCluster(n_cores)
  doParallel::registerDoParallel(cl)
  scores <- foreach(len = range) %dopar% {
    signature = all_genes[1:len]
    datasetCut <- dataset[signature, ]
    distMatrix <- Dist(t(datasetCut), nbproc=6, method="kendall")
    Hobotnica(distMatrix, annotation)
  }
  names(scores) <- range
  parallel::stopCluster(cl)
  
  return (scores)
  
}



de_deseq2_sorted <- de_deseq2[order(de_deseq2$padj),]
de_deseq2_sorted_filtered <- de_deseq2_sorted[]
scores_deseq2 <- generateScoresByLength(cm, anno$Samples,rownames(de_deseq2_sorted) )




for (score in unlist(scores_deseq2)) {
  cat(score)
  cat("\n")
}


