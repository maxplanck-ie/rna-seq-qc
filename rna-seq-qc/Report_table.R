library(gdata)

.Library
.Library.site
.libPaths()

sessionInfo()

args = commandArgs(TRUE)

## For debugging only!!! #######################################################
# setwd("/data/processing/kilpert/test/rna-seq-qc/Katarzyna/outdir_1M/project_report/")
# args = c('/data/boehm/sikora/nusser/A387/reads',
#         '/data/processing/kilpert/test/rna-seq-qc/Katarzyna/outdir_1M')
################################################################################

main_indir = args[1]
main_outdir = args[2]


## binaries ####################################################################
hostname = Sys.info()["nodename"]
if (hostname == "pc305.immunbio.mpg.de") {
  samtools_dir = "/home/kilpert/Software/samtools/samtools-1.2/"
} else {
  samtools_dir = "/package/samtools-1.2/bin/"
}


print("Generating project report...")


## get R1 read if paired ####################################################

unpair <- function(allfiles){
  files = c()
  names = c()
  for (file in allfiles) {
    name = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file))
    
    if (!name %in% names) {
      names = c(names, name)
      files = c(files, file)
    } 
  }  
  return(files)
}


## check if reads are paired ###################################################
is_paired = function(allfiles){
  names = c()
  paired = F
  for (file in allfiles) {
    name = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file))
    
    if (!name %in% names) {
      names = c(names, name)
    } else {
      paired = T
    }
  }  
  return(paired)
}



## FASTQ: total ################################################################

## If there is downsampling
if ( dir.exists(file.path(main_outdir,"FASTQ_downsampling")) ) {
  indir = file.path(main_outdir,"FASTQ_downsampling")
} else {
  indir = main_indir
}
  
if (file.exists(indir)){
  allfiles <- sort(list.files(indir, pattern="*.fastq.gz$", full.names=T))
  paired = is_paired(allfiles)
  
  num_reads = c()
  files = c()
  names = c()
  
  for (file in unpair(allfiles)) {
    name = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file))
    names = c(names, name)
    num = as.numeric( system(sprintf("zcat %s | wc -l | awk '{print $1/4}'", file), intern=T) )
    num_reads = c( num_reads, num )  
    cat(name, num, "\n")
  }
  report = data.frame(row.names=names, FASTQ=num_reads)
  report$FASTQ_perc = 100.0
} else {
  cat( paste("Error! Directory does NOT exist:\n",indir) )
  quit(save = "no", status = 1, runLast = FALSE)   # Exit 1
}

report


## trimmed FASTQ: trimmed ######################################################

indir = file.path(main_outdir,"Trim_Galore")
if ( file.exists( indir ) ){
  files <- sort(list.files(indir, pattern="*.fastq.gz$", full.names=T))
  files
  names = c()
  counts = c()
  for (file in unpair(files)) {
    name = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file))
    names = c( names, name )
    count = as.numeric( system(sprintf("zcat %s | wc -l | awk '{print $1/4}'", file), intern=T) )
    counts = c( counts, count )
    cat(paste(name, num, "\n"), sep=" ")
  }
  df_counts = data.frame(names=names, TrimGalore=counts)
  
  report = merge(report, df_counts, by.x=c("row.names"), by.y=("names"))
  rownames(report) = report$Row.names
  report$Row.names = NULL
  
  report$TrimGalore_perc = apply(report, 1, function(x) x["TrimGalore"]/x["FASTQ"]*100)
  report
} 


## BAM: mapped #################################################################
## HISAT2
indir = file.path(main_outdir,"HISAT2")
if ( file.exists( indir ) ){
  files <- sort(list.files(indir, pattern="*.bam$", full.names=T))
  ## files = sample(files) ## debugging only!!!
  files
  counts = c()
  names = c()
  
  for (file in files) {
    name = gsub(".bam$","",basename(file))
    names = c(names, name)

    ## count mapped (-F4) and primary (-F256) fragments
    count = as.numeric( system(sprintf("%s view -F260 %s | cut -f1 | sort -u | wc -l", file.path(samtools_dir,"samtools"), file), intern=T) )

    counts = c( counts, count )
    cat(paste(name, num, "\n"), sep=" ")
  }
  df_counts = data.frame(names=names, HISAT2=counts)
  
  report = merge(report, df_counts, by.x=c("row.names"), by.y=("names"))
  rownames(report) = report$Row.names
  report$Row.names = NULL
  
  report$HISAT2_perc = apply(report, 1, function(x) x["HISAT2"]/x["FASTQ"]*100)
  report
}


## TopHat2
indir = file.path(main_outdir,"TopHat2")
if ( file.exists( indir ) ){
  files <- sort(list.files(indir, pattern="*.bam$", full.names=T))
  ## files = sample(files) ## debugging only!!!
  files
  counts = c()
  names = c()
  
  for (file in files) {
    name = gsub(".bam$","",basename(file))
    names = c(names, name)
    
    if (paired==T) {
      count = as.numeric( system(sprintf("%s view -c -f3 -F260 %s | awk '{print $1/2}'", file.path(samtools_dir,"samtools"), file), intern=T) )
    } else {
      count = as.numeric( system(sprintf("%s view -c -F260 %s", file.path(samtools_dir,"samtools"), file), intern=T) )
    }
    counts = c( counts, count )
    cat(paste(name, num, "\n"), sep=" ")
  }
  df_counts = data.frame(names=names, TopHat2=counts)
  
  report = merge(report, df_counts, by.x=c("row.names"), by.y=("names"))
  rownames(report) = report$Row.names
  report$Row.names = NULL
  
  report$TopHat2_perc = apply(report, 1, function(x) x["TopHat2"]/x["FASTQ"]*100)
  report
}


## counts ######################################################################
##featureCounts
indir = file.path(main_outdir,"featureCounts")
infile = file.path(indir,"counts.txt")
if ( file.exists( infile ) ) {
 cat(infile, "\n") 
 counts = read.table(infile, header=TRUE)
 head(counts)

 df_counts = data.frame("featureCounts"=apply(counts, 2, sum))
 report = merge(report, df_counts, by.x=c("row.names"), by.y=("row.names"))
 row.names(report) = report$Row.names
 report$Row.names = NULL
 
 report$featureCounts_perc = apply(report, 1, function(x) x["featureCounts"]/x["FASTQ"]*100)
 report
}

##htseq-count
indir = file.path(main_outdir,"htseq-count")
infile = file.path(indir,"counts.txt")
if ( file.exists( infile ) ) {
  cat(infile, "\n") 
  counts = read.table(infile, header=TRUE)
  head(counts)
  
  df_counts = data.frame( "htseq_count"=apply(counts, 2, sum) )
  report = merge(report, df_counts, by.x=c("row.names"), by.y=("row.names"))
  row.names(report) = report$Row.names
  report$Row.names = NULL

  report$htseq_counts_perc = apply(report, 1, function(x) x["htseq_count"]/x["FASTQ"]*100)
  report
}

write.table(report,"Report.tsv", sep="\t", quote=FALSE, col.names=NA) # save to file

