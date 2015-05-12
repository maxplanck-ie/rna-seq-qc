## Usage: cat DESeq2.R | /package/R-3.1.0/bin/R --vanilla --quiet --args setup.tsv counts.txt 0.05 [BioMart.tsv]

library("DESeq2")

args = commandArgs(TRUE)
## Debug only! ################################################################
# setwd("/data/jenuwein/group/kilpert/140731_MeRIP_Ausma/20_DEseq/")
# args = c('/data/jenuwein/group/kilpert/140731_MeRIP_Ausma/sampleInfo.tsv',
#             '/data/jenuwein/group/kilpert/140731_MeRIP_Ausma/18_rna-seq-qc/featureCounts/counts.txt',
#             '0.05',
#             "/home/kilpert/git/rna-seq-qc/rna-seq-qc/mm10.gene_names")
###############################################################################

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

## extract only the columns specified in the sampleInfo
countdata = countdata[,sampleInfo[,1]]
head(countdata)

## check if columns names (from the count table) and sample names (from the sampleInfo) match
if ( length(setdiff(colnames(countdata), as.character(sampleInfo[,1]))) > 0 ) {
  cat("Error! Sample names in setup table and count table do NOT match!\n")
  cat(paste(sampleInfoFilePath, "\n"))
  cat(paste(paste(sampleInfo[,1], sep=" ")),"\n")
  cat(paste(countFilePath, "\n"))
  cat(paste(paste(colnames(countdata), sep=" ")),"\n")
  q("no", 1, FALSE) # exit code 1
}

# check if sample names are the same in the input files
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
summary(res)

################################################################################
## gene names dict if available
################################################################################

if (file.exists(geneNamesFilePath)) { 
  cat(paste("Gene names file found\n")) 
  geneNames = read.csv(geneNamesFilePath, sep="\t", header=F, row.names=1, stringsAsFactors=FALSE)
  
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
write.table(de_total[order(de_total$padj, decreasing=F),],"DESeq2.all.tsv", sep="\t", quote=FALSE, col.names=NA)

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

# topN genes by pvalue
d = data.frame(id=rownames(res), padj=res$padj)
d_topx_padj = d[order(d$padj, decreasing=F),][1:topN,]
d_topx_padj
plotdata = assay(rld)[d_topx_padj$id,]

if ( exists("gene_names_dic") ) rownames(plotdata) = id_to_gene_name(rownames(plotdata))  # exchange ids by gene names

pdf(sprintf("Fig6.gene_clustering_top%i_DE_genes.pdf",topN), pointsize = 9)
heatmap.2(plotdata, scale="row", trace="none", dendrogram="column",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
          main=sprintf("Top %d DE genes (by p-value)", topN), keysize=1,
          margins = c(10,12))
dev.off()


################################################################################

# report on versions used
sink("DESeq2.session_info.txt")
sessionInfo()
sink()

