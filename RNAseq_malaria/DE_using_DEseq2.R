library(DESeq2)
# input the raw count data
directory <- "../2_Read_counts/"
sampleFiles <- list.files(directory)
sampleCondition <- sub("_.{2,}","",sampleFiles)
sampleNames <- sub("_\\d.{1,}count","",sampleFiles) # / has to escape itself, when using it to escape character.
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)
summary(sampleTable)

ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition )

#filtering the low count genes
keep <- rowSums(counts(ddsHTseq)) >= 10
ddsHTseq <- ddsHTseq[keep,]
ddsHTseq$condition

# Differential expression analysis
dds <- DESeq(ddsHTseq)
res <- results(dds, contrast = c("condition","D8","D20"))
res
summary(res) #dow and up regulated genes are too few, too much low counts.
##the meaning of col name
mcols(res)

#Log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "condition_D8_vs_D20",res = res)


#custom res table
#p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),]
resOrdered

# filtering sum
sum(res$padj < 0.05, na.rm = TRUE) # P value of FDR 
sum(res$pvalue < 0.05, na.rm = TRUE)
# FDR cutoff
res05 <- results(dds, alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm = TRUE) # why those two values are not the same?


# MA-plot 
# with out stastical moderation
plotMA(res, ylim=c(-5,5)) #the plot means almost no change
plotMA(resLFC, ylim=c(-2,2)) #  remove the noise associated with log2 fold  changes from low count genes without requiring arbitrary filtering threshold


resSig <- subset(res, padj < 0.05)
summary(resSig)
write.table(resSig, file = "../3_DE/diffExp.0.05.tab", sep= "\t", quote = FALSE)

#Quality assessment
#Expression heatmap 
library("pheatmap")

ntd <- normTransform(ddsHTseq)
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds,blind = T)
select <- order(rowMeans(counts(dds, normalized= TRUE)), decreasing = TRUE)
df <- as.data.frame(colData(dds)["condition"])
pheatmap(assay(rld)[select,],cluster_rows=FALSE,cluster_cols =FALSE,,annotation_col=df)

pheatmap(resSig)
pheatmap(ddsHTseq)
pdf("../3_DE/PCA.pdf")
plotPCA(rld, intgroup=c("condition"))
dev.off()

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- with(colData(dds), paste(vsd$condition, sampleFiles, sep = ":"))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

