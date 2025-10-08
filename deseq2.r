args <- commandArgs(trailingOnly = TRUE)

require(DESeq2)

countData <- read.csv(args[1], header = TRUE, sep = ",")
metaData <- read.csv(args[2], header = TRUE, sep = ",")
dds <- DESeqDataSetFromMatrix(countData=countData, colData=metaData, design=~y, tidy=T)
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, args[3])