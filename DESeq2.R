library(DESeq2)

cts1 <- as.matrix(read.csv("data/SRR11497577_featureCounts.txt", sep='\t', header=TRUE, row.names="Geneid", comment="#"))
cts1 = cts1[,colnames(cts1) == "HISAT2.SRR11497577.bam", drop=FALSE]
colnames(cts1) <- "SRR11497577"

cts2 <- as.matrix(read.csv("data/SRR11497582_featureCounts.txt", sep='\t', header=TRUE, row.names="Geneid", comment="#"))
cts2 = cts2[,colnames(cts2) == "HISAT2.SRR11497582.bam", drop=FALSE]
colnames(cts2) <- "SRR11497582"

cts3 <- as.matrix(read.csv("data/SRR11497581_featureCounts.txt", sep='\t', header=TRUE, row.names="Geneid", comment="#"))
cts3 = cts3[,colnames(cts3) == "HISAT2.SRR11497581.bam", drop=FALSE]
colnames(cts3) <- "SRR11497581"

cts4 <- as.matrix(read.csv("data/SRR11497578_featureCounts.txt", sep='\t', header=TRUE, row.names="Geneid", comment="#"))
cts4 = cts4[,colnames(cts4) == "HISAT2.SRR11497578.bam", drop=FALSE]
colnames(cts4) <- "SRR11497578"

cts <- cbind(cts1, cts2, cts3, cts4)

colData <- data.frame(specialization = c("naive", "central_memory", "naive", "central_memory"), row.names = c("SRR11497577", "SRR11497582", "SRR11497581", "SRR11497578"))

colData$specialization <- as.factor(colData$specialization)

all(colnames(cts) %in% rownames(colData))

all(colnames(cts) == rownames(colData))

numeric_cts <- apply(cts, 2, as.numeric)

rownames(numeric_cts) <- rownames(cts)

filtered_cts <- numeric_cts[rowSums(numeric_cts) > 10, ]

dds <- DESeqDataSetFromMatrix(countData = filtered_cts,
                              colData = colData,
                              design = ~ specialization)

dds

dds$specialization <- relevel(dds$specialization, ref = "naive")

dds <- DESeq(dds)

res <- results(dds)

res

summary(res)

res0.01 <- results(dds, alpha = 0.01)

resultsNames(dds)

# results(dds, constrast = c("specialization", "central_memory", "effector_memory"))

plotMA(res0.01)

res0.01 <- res0.01[order(res0.01$padj),]

par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000077984", intgroup="specialization")
plotCounts(dds, gene="ENSG00000113088", intgroup="specialization")
plotCounts(dds, gene="ENSG00000181847", intgroup="specialization")
plotCounts(dds, gene="ENSG00000116741", intgroup="specialization")
plotCounts(dds, gene="ENSG00000134193", intgroup="specialization")
plotCounts(dds, gene="ENSG00000132965", intgroup="specialization")

geneList=c(which(results(dds, tidy=TRUE, alpha=0.01)[,7]<=2.e-09))
n=length(geneList)
geneNames=attributes(res[1:n,])$rownames
