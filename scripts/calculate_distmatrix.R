library("amap")
args = commandArgs(trailingOnly=TRUE)


countMatrixFile <- args[1]

signatureFile <- args[2]

signature <- readLines(signatureFile)

#print(signature)


cm <- read.table(countMatrixFile, header=T, sep=",")
cm <- data.frame(cm[, -1], row.names=cm[, 1])

cm_subset <- cm[signature, ]
cat("SUBSET ShAPE: ")
cat(dim(cm_subset))
cat("\n")
distMatrix <- Dist(t(cm_subset), method="kendall", nbproc=6)
write.csv(as.data.frame(as.matrix(distMatrix)), file=paste0(countMatrixFile, ".distmatrix"))
