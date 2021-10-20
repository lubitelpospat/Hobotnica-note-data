
library(amap)
library(edgeR)
library(biomaRt)
library(Hobotnica)
library(MASS)
args = commandArgs(trailingOnly=TRUE)

countMatrixFile <- args[1]
annotationFile <- args[2]
mode <- args[3]



signature <- c('ACTR3B', 'ANLN', 'BAG1', 'BCL2', 'BIRC5', 'BLVRA', 'CCNB1', 'CCNE1', 'CDC20', 'CDC6', 'CDH3', 'CENPF', 'CEP55', 'CXXC5', 'EGFR', 'ERBB2', 'ESR1', 'EXO1', 'FGFR4', 'FOXA1', 'FOXC1', 'GPR160', 'GRB7', 'KIF2C', 'KRT14', 'KRT17', 'KRT5', 'MAPT', 'MDM2', 'MELK', 'MIA', 'MKI67', 'MLPH', 'MMP11', 'MYBL2', 'MYC', 'NAT1', 'NDC80', 'NUF2', 'ORC6L', 'PGR', 'PHGDH', 'PTTG1', 'RRM2', 'SFRP1', 'SLC39A6', 'TMEM45B', 'TYMS', 'UBE2C', 'UBE2T')


print("SIGNATURE LENGTH")

print(length(signature))
if (mode == "ensembl") {
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
signature_ensembl <- getBM(values=signature,
			    filters = "external_gene_name",
			    mart = mart,
  attributes = c("external_gene_name", "entrezgene_id",
    "hgnc_symbol", "description",
    "chromosome_name", "strand", "ensembl_gene_id"))$ensembl_gene_id

print("ENSEMBL PAM50 Signature")
} else if (mode != "gid") {
	stop ("mode should be one of: `ensembl`, `gid`. Exiting")
}

#print(signature_ensembl)


cm <- read.table(countMatrixFile, header=T, sep=",")
cm <- data.frame(cm[, -1], row.names=cm[, 1])
print("COUNT MATRIX SHAPE")
print(dim(cm))
anno_raw = read.table(annotationFile, header=TRUE, sep=",")

colnames(anno_raw) <- c("Run", "source_name")
anno = data.frame(anno_raw$source_name, row.names = anno_raw$Run)
colnames(anno) <- "group"
print(head(anno))


#print(colnames(cm))
cat("SHARED SAMPLES SIZE: ")
cat(length(intersect(colnames(cm), rownames(anno))))
cat("\n")

#cm <- cm[, rownames(anno)]
#cm_cpm <- cpm(cm)
cm_cpm <- cm #cpm(cm)
randomSignaturesList <- list()
randomSigScores <- list()
for (i in 1:3000) {
            randomSignaturesList[[paste0("random_", i)]] <- sample(rownames(cm_cpm), length(signature))
}
for (name in names(randomSignaturesList)) {
                set_filt <- randomSignaturesList[[name]]
                xxx <- cm_cpm[set_filt, ]
                distMatrix <- Dist(t(xxx), method="kendall", nbproc=10)
                
                randomSigScores[[name]] <- Hobotnica(distMatrix, anno$group)
}
	

if (mode == "ensembl") {
	set_filt <- intersect(rownames(cm_cpm), signature_ensembl)
} else  if (mode == "gid"){
	set_filt <- intersect(rownames(cm_cpm), signature)
} else {
	stop("FUBAR")
}
cm_cpm <- cm_cpm[set_filt, ]
print("CM SUBSET SHAPE")
print(dim(cm_cpm))
distMatrix <- as.matrix(Dist(t(cm_cpm), method="kendall", nbproc=10))
score <- Hobotnica(distMatrix, anno$group)
write.csv(as.data.frame(as.matrix(distMatrix)), file=paste0(countMatrixFile, ".pam50.distmatrix"))
cat("SCORE: ")
cat(score)
cat("\n")
pval <- min(1, (length(which(unlist(randomSigScores) >= score)) + 1)/length(unlist(randomSigScores)))
cat("P-value: ")
cat(pval)
cat("\n")
