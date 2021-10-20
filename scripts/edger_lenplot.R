library("edgeR")
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


y <- DGEList(counts=cm,group=as.factor(anno$Samples))
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
cor <- topTags(et, n = dim(et)[1], sort.by='none')
de_edger <- cor$table

de_edger_sorted <- de_edger[order(de_edger$FDR),]




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



scores_edger <- generateScoresByLength(cm, anno$Samples,rownames(de_edger_sorted) )




for (score in unlist(scores_edger)) {
  cat(score)
  cat("\n")
}


