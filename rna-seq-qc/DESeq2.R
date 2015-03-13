## Usage: cat DESeq2.R | R --vanilla --args count_folder sampleInfo.csv 0.05 [BioMart.tsv]

library("DESeq2")

args = commandArgs(TRUE)

print("Running DESeq2 from rna-seq-qc...")

## FDR significance threshold
fdr = as.numeric(args[3])
if ( is.na(fdr) ) fdr = 0.05  # default FDR

#from command line
sampleInfoFilePath = args[1]
countFilePath = args[2]
biomartFilePath = args[4]   # BioMart file with ensembl and symbol names

cat(paste("Sample info CSV:", sampleInfoFilePath, "\n"))
cat(paste("Count file:", countFilePath, "\n"))
cat(paste("FDR:", fdr, "\n"))
cat(paste("BioMart file:", biomartFilePath, "\n"))

## sampleInfo (setupt of the experiment)
sampleInfo = read.table(sampleInfoFilePath, header=TRUE)
sampleInfo = DataFrame(sampleInfo)
sampleInfo

## count matrix (e.g. from DESeq or featureCounts)
countdata = read.table(countFilePath, header=TRUE)
countdata = DataFrame(countdata)
head(countdata)

dds = DESeqDataSetFromMatrix(
  countData = countdata,
  colData = sampleInfo,
  design = ~ condition)
dds
colnames(dds) = sampleInfo$name
head(assay(dds))

################################################################################
## Extra data collecting some measures for every sample (e.g. total
## counts, scaling factors,...)
info <- data.frame(row.names=sampleInfo$name)
################################################################################
## counts per sample
apply(assay(dds), 2, sum) 
info$total_counts = apply(assay(dds), 2, sum)   # add to info
info

## count table
head(assay(dds))
write.table(assay(dds),"DESeq2.counts.tsv", sep="\t", quote=FALSE) # save to file

## DE analysis
assign("last.warning", NULL, envir = baseenv())
dds = DESeq(dds)
warnings()
sink("DESeq2.WARNINGS.txt"); warnings(); sink() # save warnings to file

## show size factors used for read count normalisation
sizeFactors(dds)
info$size_factors = sizeFactors(dds) 
info

## dispersion plot
pdf("Fig1.dispersion_plot.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()

## get results
res = results(dds)
head(res)

################################################################################
## extend results by gene names obtained via Ensembl IDs from BioMart
################################################################################
if (file.exists(biomartFilePath)) { 
  cat(paste("BioMart.tsv found\n")) 
  bmGeneNames = read.csv(biomartFilePath, sep="\t", header=TRUE, row.names=1)
  
  res_output = merge(res,
                        bmGeneNames,
                        by.x = "row.names",
                        by.y = "ensembl_gene_id",
                        all.x = TRUE) # set gene name to "NA" if not available
  rownames(res_output) = res_output$Row.names
  res_output$Row.names = NULL
  res@listData = c(res@listData, external_gene_id=list(res_output$external_gene_id))  # append new list to res@listData; keep everything in one place
  head(res_output)
  tail(res_output) 
} 
else {
    cat(paste("BioMart.tsv NOT found\n"))
#   ## download if necessary
#   library("biomaRt")
#   cat(paste("BioMart.tsv NOT found. Downloading...\n")) 
#   ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#   bmGeneNames = getBM(attributes = c("ensembl_gene_id", "external_gene_id"),
#                       filters = "ensembl_gene_id",
#                       values = rownames(res),
#                       mart = ensembl)
#   write.table(bmGeneNames,"BioMart.tsv", sep="\t", quote=FALSE, col.names=NA) 
}

################################################################################

## DE
##fdr = 0.05
de_total = res[which(res$padj < fdr),]
length(de_total[,1])

de_up = de_total[which(de_total$log2FoldChange>0),]
de_up = de_up[order(de_up$log2FoldChange, decreasing=T),]   # order by log2FoldChange
length(de_up[,1])
write.table(de_up,"DESeq2.up.tsv", sep="\t", quote=FALSE, col.names=NA)

de_down = de_total[which(de_total$log2FoldChange<0),]
de_down = de_down[order(de_down$log2FoldChange),]           # order by log2FoldChange
length(de_down[,1])
write.table(de_down,"DESeq2.down.tsv", sep="\t", quote=FALSE, col.names=NA)

# save info to metrics file
write.table(info,"DESeq2.metrics.tsv", sep="\t", quote=FALSE, col.names=NA)
 
# MA plot
pdf("Fig2.MA-plot.pdf")
plotMA(res, alpha=fdr, ylim=c(-2,2), 
       main=sprintf("MA-plot\n(FDR: %.2f, up: %d, down: %d)",fdr,length(de_up[,1]),length(de_down[,1])),
       ylab="log2 fold change")
dev.off()

## Histogram of p-values
pdf("Fig3.p-values_histogram.pdf")
hist(res$pvalue, breaks=20, col="grey", main="Histogram of p-values", xlab="p-value")
dev.off()

## Histogram of adjusted p-values
pdf("Fig4.padj_histogram.pdf")
hist(res$padj, breaks=20, col="grey", main="Histogram of adjusted p-values", xlab="padj")
abline(v=fdr, col="red", lwd=1)
dev.off()

# ## Independent filtering
attr(res,"filterThreshold")
plot(attr(res,"filterNumRej"), type="b",
     ylab="number of rejections",
     xlab="quantiles of mean of normalized counts")
 
plot(res$baseMean+1, -log10(res$padj),
     log="x", 
     xlab="mean of normalized counts",
     cex=.4, col=rgb(0,0,0,.3))
abline(h=-log10(fdr), col="red", lwd=1)

################################################################################
## rlog transform; for clustering and ordination (e.g PCA)
rld = rlog(dds)
head(assay(rld))
# show that DEseq's rlog works; not really needed for data analysis
##par(mfrow=c(1,2))
##plot(log2(1+counts(dds, normalized=T)[,1:2]), col="#00000020", pch=20, cex=0.3)     #log2
##plot(assay(rld)[,1:2], col="#00000020", pch=20, cex=0.3)                            #rlog (DESeq); is superior

## Sample distances
sampleDists = dist(t(assay(rld)))
sampleDists

## Euclidean sample distance heatmap
sampleDistMatrix = as.matrix(sampleDists)
sampleDistMatrix
rownames(sampleDistMatrix) = sprintf("%s %s", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
colnames(sampleDistMatrix) = sprintf("%s %s", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
sampleDistMatrix
library("gplots")
library("RColorBrewer")
colours = colorRampPalette(rev(brewer.pal(9, "GnBu")))(255)
pdf("Fig5.Heatmap.pdf")
heatmap.2(sampleDistMatrix,trace="none",col=colours,
          main="Heatmap (Euclidean distances)",
          keysize=1.3, margins=c(8,8))
dev.off()

## PCA
pdf("PCA.pdf")
print(plotPCA(rld, intgroup=c("condition")))
dev.off()

## Gene clustering (across samples)
library("genefilter")
n=50
topVarGenes = head(order(rowVars(assay(rld)), decreasing=T), n)
pdf(sprintf("Fig7.gene_clustering_top%i.pdf",n))
heatmap.2(assay(rld)[topVarGenes, ], scale="row", trace="none", dendrogram="column",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
          main=sprintf("Gene clustering (top %d)", n),keysize=1.3,
          margins = c(4,12))
dev.off()
# ################################################################################

# report on versions used
sink("DESeq2.session_info.txt")
sessionInfo()
sink()


