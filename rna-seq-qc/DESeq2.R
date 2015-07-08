## Usage: cat DESeq2.R | /package/R-3.2.0/bin/R --vanilla --quiet --args setup.tsv counts.txt 0.05 species.gene_names

.Library
.Library.site
.libPaths()

library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")

sessionInfo()

args = commandArgs(TRUE)

# ## For debugging only!!! #######################################################
# setwd("/data/processing/kilpert/test/rna-seq-qc/Ausma/PE_mm10_subset_hisat/DESeq2/")
# args = c('/data/manke/kilpert/datasets/Ausma/subset/sampleInfo.tsv',
#          '/data/processing/kilpert/test/rna-seq-qc/Ausma/PE_mm10_subset_hisat/featureCounts/counts.txt',
#          '0.05',
#          '/home/kilpert/git/rna-seq-qc/rna-seq-qc/mm10.gene_names')
# ################################################################################

plotVolcano <- function(res_obj, data=plot) {
  # Volcano plot
  xlim = c(-2.5,2.5)
  ylim = c(0,20)
  cex=c(0.3,0.5)
  plotdata = data.frame(log2FoldChange=res_obj$log2FoldChange, padj=res_obj$padj )
  plotdata = plotdata[!is.na(plotdata),]
  plotdata$cex = cex[[1]]
  plotdata$pch = 19
  plotdata$col = "#525252"
  plotdata$col[plotdata$padj<=fdr] = "#cd0000"
  
  plotdata$pch[plotdata$log2FoldChange<xlim[[1]]] = 5
  plotdata$cex[plotdata$log2FoldChange<xlim[[1]]] = cex[[2]]
  plotdata$log2FoldChange[plotdata$log2FoldChange<xlim[[1]]] = xlim[[1]]

  plotdata$pch[plotdata$log2FoldChange>xlim[[2]]] = 5
  plotdata$cex[plotdata$log2FoldChange>xlim[[2]]] = cex[[2]]
  plotdata$log2FoldChange[plotdata$log2FoldChange>xlim[[2]]] = xlim[[2]]
  
  plotdata$pch[-log10(plotdata$padj) > ylim[[2]]] = 2
  plotdata$cex[-log10(plotdata$padj) > ylim[[2]]] = cex[[2]]
  plotdata$padj[-log10(plotdata$padj) > ylim[[2]]] = 10^-ylim[[2]]
  
  #head(plotdata)
  #dim(plotdata)
  plot(plotdata$log2FoldChange, -log10(plotdata$padj),
       main=sprintf("Volcano plot\n(FDR: %.2f, up: %d, down: %d)",fdr,length(de_up[,1]),length(de_down[,1])),
       xlab="log2-fold change",
       ylab="-log10 q-value",
       xlim=xlim,
       ylim=ylim,
       cex=plotdata$cex, pch=plotdata$pch,
       col=plotdata$col)
  abline(h=-log10(fdr), col=alpha(4,0.5), lwd=4)
  abline(v=0, col=alpha(2,0.5), lwd=4)
}

################################################################################


print("Running DESeq2 from rna-seq-qc...")

## FDR significance threshold
fdr = as.numeric(args[3])
if ( is.na(fdr) ) fdr = 0.05  # default FDR

topN =as.numeric(args[5])
if ( is.na(topN) ) topN = 50  # use topN genes for plots

#from command line
sampleInfoFilePath = args[1]
countFilePath = args[2]
geneNamesFilePath = args[4]   # BioMart file with ensembl and symbol names

cat(paste("Sample info CSV:", sampleInfoFilePath, "\n"))
cat(paste("Count file:", countFilePath, "\n"))
cat(paste("FDR:", fdr, "\n"))
cat(paste("Gene names:", geneNamesFilePath, "\n"))
cat(paste("Number of top N genes:", topN, "\n"))

## sampleInfo (setupt of the experiment)
sampleInfo = read.table(sampleInfoFilePath, header=TRUE)
sampleInfo = DataFrame(sampleInfo)
sampleInfo = sampleInfo[order(sampleInfo$name, decreasing=F),]  # order by sample name
sampleInfo

## count matrix (e.g. from DESeq or featureCounts)
countdata = read.table(countFilePath, header=TRUE)
countdata = DataFrame(countdata)
countdata = countdata[,order(names(countdata), decreasing=F)]  # order column names
head(countdata)

## 1st check: if names of the setup table are subset of the count matrix column names
if ( ! all( is.element(sampleInfo[,1], colnames(countdata)) ) ) {
  cat("Error! Count table column names and setup table names do NOT match!\n")
  print(as.character(sampleInfo[,1]))
  print(colnames(countdata))
  quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
}

## extract only the columns specified in the setup table
countdata = countdata[,as.character(sampleInfo[,1])]
head(countdata)

## 2nd: check if the ORDER of sample names matches as well
if ( ! all(as.character(sampleInfo$name) == colnames(countdata)) ) {
  cat("Error! Count table column names and setup table names do NOT match!\n")
  print(as.character(sampleInfo$name))
  print(colnames(countdata))
  quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
}

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
write.table(assay(dds),"DESeq2.counts_raw.tsv", sep="\t", quote=FALSE, col.names=NA) # save to file

## DE analysis
assign("last.warning", NULL, envir = baseenv())
dds = DESeq(dds)
warnings()
sink("DESeq2.WARNINGS.txt"); warnings(); sink() # save warnings to file

## show size factors used for read count normalisation
sizeFactors(dds)
info$size_factors = sizeFactors(dds) 
info

# save normalized counts to file
write.table(counts(dds, normalized=T),"DESeq2.counts_normalized.tsv", sep="\t", quote=FALSE, col.names=NA)

# ## size factors plot
# plotdata = data.frame(name=colnames(dds), size_factor=sizeFactors(dds))
# plotdata
# ggplot(plotdata, aes(x=name)) +
#   geom_bar(aes(weight=size_factor), fill="darkseagreen", width=.7) +
#   theme_bw(base_size=14) +
#   geom_text(aes(y=size_factor, label=sprintf("%.3f",size_factor)), size = 4, hjust = 0.5, vjust = -1) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("") +
#   ylab("Size factor")
# ggsave(file=sprintf("Fig7.Size_factors.pdf"))

## dispersion plot
pdf("Fig1.dispersion_plot.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()

## get results
res = results(dds)
head(res)
summary(res)
dim(res)

################################################################################
## gene names dict if available
################################################################################

if (file.exists(geneNamesFilePath)) { 
  cat(paste("Gene names file found\n")) 
  #geneNames = read.csv(geneNamesFilePath, sep="\t", header=F, row.names=1, stringsAsFactors=FALSE)
  geneNames = read.csv(geneNamesFilePath, sep="\t", header=F, stringsAsFactors=FALSE)
  geneNames = geneNames[!duplicated(geneNames[,1]),]
  rownames(geneNames) = geneNames[,1]
  geneNames[,1] = NULL
  head(geneNames)
  
  if (length( intersect( gsub("\\..*", "", res@rownames), rownames(geneNames) ) ) > 0) {
    cat(paste("Names matching to IDs found\n")) 
    
    ## make a dictionary
    gene_names_dic = geneNames[[1]]
    names(gene_names_dic) = rownames(geneNames)
    ##gene_names_dic["ENSMUSG00000025332"]
  }
}

## generate a dataframe from ids
id_to_gene_name = function(ids) {
  d = data.frame(IDs=gsub("\\..*", "", ids), gene_names=NA)
  d$gene_names = gene_names_dic[ as.character(d$IDs) ]
  head(d)
  
  # some might be NAs; replace those by original ID
  d[which(is.na(d$gene_names)),]$gene_names = as.character(d[which(is.na(d$gene_names)),]$IDs)
  head(d)
  return(d$gene_name)
}

################################################################################

if ( exists("gene_names_dic") ) {
  cat("Gene names are available\n")
  # update res with gene names
  res@listData$gene_names = id_to_gene_name(rownames(res))
}


## DE ##########################################################################
de_total = res[which(res$padj < fdr),]
length(de_total[,1])
write.table(de_total[order(de_total$padj, decreasing=F),],"DESeq2.de_all.tsv", sep="\t", quote=FALSE, col.names=NA)

de_up = de_total[which(de_total$log2FoldChange>0),]
de_up = de_up[order(de_up$padj, decreasing=F),]   # order by adjusted p-value
length(de_up[,1])
write.table(de_up,"DESeq2.de_up.tsv", sep="\t", quote=FALSE, col.names=NA)

de_down = de_total[which(de_total$log2FoldChange<0),]
de_down = de_down[order(de_down$padj, decreasing=F),]           # order by adjusted p-value
length(de_down[,1])
write.table(de_down,"DESeq2.de_down.tsv", sep="\t", quote=FALSE, col.names=NA)

# save info to stats file
write.table(info,"DESeq2.stats.tsv", sep="\t", quote=FALSE, col.names=NA)

# MA and volcano plot
pdf("Fig2.MA_and_Volcano_plot.pdf", width=12, height=6)
par(mfrow=c(1,2))
plotMA(res, alpha=0.1, ylim=c(-2,2), 
       main=sprintf("MA-plot\n(FDR: %.2f, up: %d, down: %d)",fdr,length(de_up[,1]),length(de_down[,1])),
       ylab="log2 fold change")
plotVolcano(res)
dev.off()

# ## Histogram of p-values
# pdf("FigX.p-values_histogram.pdf")
# hist(res$pvalue, breaks=20, col="grey", main="Histogram of p-values", xlab="p-value")
# dev.off()

## Histogram of adjusted p-values
pdf("Fig3.padj_histogram.pdf")
hist(res$padj, breaks=20, col="grey", main="Histogram of adjusted p-values", xlab="padj")
abline(v=fdr, col="red", lwd=1)
dev.off()

# # ## Independent filtering
# attr(res,"filterThreshold")
# plot(attr(res,"filterNumRej"), type="b",
#      ylab="number of rejections",
#      xlab="quantiles of mean of normalized counts")
#  
# plot(res$baseMean+1, -log10(res$padj),
#      log="x", 
#      xlab="mean of normalized counts",
#      ylab="-log10 padj",
#      cex=.4, col=rgb(0,0,0,.3))
# abline(h=-log10(fdr), col="red", lwd=1)

################################################################################
## rlog transform; for clustering and ordination (e.g PCA)
rld = rlog(dds)
head(assay(rld))

# save rlog tranformed counts to file
write.table(assay(rld),"DESeq2.counts_rlog.tsv", sep="\t", quote=FALSE, col.names=NA)

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
rownames(sampleDistMatrix) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
colnames(sampleDistMatrix) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
sampleDistMatrix

colours = colorRampPalette(rev(brewer.pal(9, "GnBu")))(255)
pdf("Fig4.Heatmap.pdf", width=6, height=6)
heatmap.2(sampleDistMatrix,trace="none",col=colours,
          main="Heatmap\n(Euclidean distances)",
          keysize=1.2,
          cex.main=3,
          cexRow=0.8, cexCol=0.8, margins=c(8,8),
          cellnote=round(sampleDistMatrix,1),
          notecol="black")
dev.off()

## PCA
data <- plotPCA(rld, intgroup=c("name", "condition"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=name, shape=condition)) +
  geom_hline(aes(yintercept=0), colour="grey") +
  geom_vline(aes(vintercept=0), colour="grey") +
  geom_point(size=5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA\n")
ggsave(file=sprintf("Fig5.PCA.pdf"), width=7, height=6)

##pdf("PCA.pdf")
##plotPCA(rld, intgroup=c("name", "condition"))
##dev.off()

# topN genes by pvalue
d = data.frame(id=rownames(de_total), padj=de_total$padj)
if ( length(rownames(d)) < topN ) topN = length(rownames(d))

d_topx_padj = d[order(d$padj, decreasing=F),][1:topN,]
d_topx_padj
plotdata = assay(rld)[as.character(d_topx_padj$id),]  # <- error
plotdata

## test
setdiff( as.character(d_topx_padj$id), rownames(plotdata))


# rownames(plotdata) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")
# colnames(plotdata) = sprintf("%s\n(%s)", colnames(rld), rld$condition) #paste(colnames(rld), rld$condition, sep="-")

if ( exists("gene_names_dic") ) rownames(plotdata) = id_to_gene_name(rownames(plotdata))  # exchange ids by gene names
plotdata

pdf(sprintf("Fig6.gene_clustering_top%i_DE_genes.pdf",topN), pointsize = 9)
heatmap.2(plotdata, scale="row", trace="none", dendrogram="column",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
          main=sprintf("Top %d DE genes (by p-value)", topN), keysize=1,
          margins = c(8,10),
          cexRow=0.7, cexCol=0.9)
dev.off()


################################################################################

# report on versions used
sink("DESeq2.session_info.txt")
sessionInfo()
sink()

