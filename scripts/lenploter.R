library("DESeq2")
library("edgeR")
library("limma")
library("foreach")
library("doParallel")
library("amap")
library("Hobotnica")
PADJ_BORDER <- 0.05
LOGFC_BORDER <- 2
args = commandArgs(trailingOnly=TRUE)

countMatrixFile <- args[1]
annotationFile <- args[2]
deMethod <- args[3]

generateScoresByLength <- function(dataset, annotation, all_genes, range=10:200, n_cores=2) {
  
  cl <- parallel::makeForkCluster(n_cores)
  doParallel::registerDoParallel(cl)
  scores <- foreach(len = range) %dopar% {
    signature = all_genes[1:len]
    set_filt <- intersect(signature, rownames(cm))
    if (length(set_filt) == 0) {
    	stop("Genes from signature and from count matrix have a zero intersection")
    }
    datasetCut <- dataset[set_filt, ]
    distMatrix <- Dist(t(datasetCut), nbproc=6, method="kendall")
    Hobotnica(distMatrix, annotation)
  }
  names(scores) <- range
  parallel::stopCluster(cl)
  
  return (scores)
  
}


cm <- read.table(countMatrixFile, header=T, sep=",")
cm <- data.frame(cm[, -1], row.names=cm[, 1])
#print("COUNT MATRIX SHAPE")
#print(dim(cm))
anno_raw = read.table(annotationFile, header=TRUE, sep=",")

colnames(anno_raw) <- c("Run", "source_name")
anno = data.frame(anno_raw$source_name, row.names = anno_raw$Run)
colnames(anno) <- "Samples"
anno$Samples <- as.factor(anno$Samples)
de_res <- data.frame()
signature <- c()
if (deMethod == "deseq2") {
	coldata <- data.frame(condition=as.factor(anno$Samples))
	dds <- DESeqDataSetFromMatrix(countData = cm,
                              colData = coldata,
                              design = ~ condition)


	cds = estimateSizeFactors(dds);
	cds = estimateDispersions(cds);
	cds = nbinomWaldTest(cds);

	de = data.frame(results(cds));
	de_deseq2 = de[!is.na(de$padj), ];
	de_res <- data.frame(padj=de_deseq2$padj, logFC=de_deseq2$log2FoldChange)
	rownames(de_res) <- rownames(de_deseq2)
} else if (deMethod=="edger") {
	y <- DGEList(counts=cm,group=as.factor(anno$Samples))
	y <- calcNormFactors(y)
	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	et <- exactTest(y)
	cor <- topTags(et, n = dim(et)[1], sort.by='none')
	de_edger <- cor$table
	de_res <- data.frame(padj=de_edger$FDR, logFC=de_edger$logFC)
	rownames(de_res) <- rownames(de_edger)
} else if (deMethod=="limma") {
		design = matrix(ncol = 2, nrow = length(anno$Samples))
		design <- as.data.frame(design)
		design[,1] <- as.numeric(rep(1,length(anno$Samples)))
		design[,2] <- as.numeric(anno$Samples)



		v <- voom(cm,design,plot=F,normalize="quantile")
		fit <- lmFit(v,design)
		fit <- eBayes(fit)
		de_limma <- topTable(fit, n = dim(fit)[1],coef = "V2", sort.by = "none", adjust="BH")
		de_res <- data.frame(padj=de_limma$adj.P.Val, logFC=de_limma$logFC)
		rownames(de_res) <- rownames(de_limma)
} else {
	stop(paste0("Invalid de method! Expected one of: deseq2, edger, limma, got ", deMethod))
}


#cat("de_res SHAPE")
#cat(dim(de_res))
#cat("\n")


lenplots <- list()
#padj_filter 

de_res_padj_filter <- de_res[de_res$padj < PADJ_BORDER, ]

lenplots[['de_res_padj_filter']] <- rownames(de_res_padj_filter)

#padj_sort 


#de_res_padj_sort <- de_res[order(de_res$padj), ]

#lenplots[['de_res_padj_sort']] <- rownames(de_res_padj_sort)

# logfc_filter 


de_res_logfc_filter <- de_res[abs(de_res$logFC) > LOGFC_BORDER,]

lenplots[['de_res_logfc_filter']] = rownames(de_res_logfc_filter)

#padj_filter_sort

de_res_padj_filter_sort <- de_res_padj_filter[order(de_res_padj_filter$padj), ]

lenplots[['de_res_padj_filter_sort']] = rownames(de_res_padj_filter_sort)

#padj_filter_padj_logfc_sort


de_res_padj_filter_padj_logfc_sort <- de_res_padj_filter[order(de_res_padj_filter$padj, -abs(de_res_padj_filter$logFC)),]

lenplots[['de_res_padj_filter_sort']] = rownames(de_res_padj_filter_sort)

#padj_filter_logfc_padj_sort 


de_res_padj_filter_logfc_padj_sort <- de_res_padj_filter[order(-abs(de_res_padj_filter$logFC), de_res_padj_filter$padj), ]

lenplots[['de_res_padj_filter_logfc_padj_sort']] = rownames(de_res_padj_filter_logfc_padj_sort)


#padj_filter_logfc_sort 



de_res_padj_filter_logfc_sort <- de_res_padj_filter[order(-abs(de_res_padj_filter$logFC)),]


lenplots[['de_res_padj_filter_logfc_sort']] = rownames(de_res_padj_filter_logfc_sort)


#logfc_filter_padj_sort

de_res_logfc_filter_padj_sort <- de_res_logfc_filter[order(de_res_logfc_filter$padj), ]
lenplots[['de_res_logfc_filter_padj_sort']] = rownames(de_res_logfc_filter_padj_sort)








for (name in names(lenplots)) {
	cat(name)
	cat("\n")
	cat("LENGTH: ")
	cat(length(lenplots[[name]]))
	cat("\n")
	#print(lenplots[[name]])
	maxLen <- min(200, length(lenplots[[name]]))
	lenplots[[name]] = generateScoresByLength(cm, anno$Samples, lenplots[[name]], range=1:maxLen, n_cores=8)
	for (score_ind in names(lenplots[[name]])) {
		cat(lenplots[[name]][[score_ind]])
		cat("\n")
	}
}











