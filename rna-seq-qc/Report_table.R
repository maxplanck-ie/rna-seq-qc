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
  samples = c()
  for (file in allfiles) {
    sample = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file))
    
    if (!sample %in% samples) {
      samples = c(samples, sample)
      files = c(files, file)
    } 
  }  
  return(files)
}


## check if reads are paired ###################################################
is_paired = function(allfiles){
  samples = c()
  paired = F
  for (file in allfiles) {
    sample = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file))
    
    if (!sample %in% samples) {
      samples = c(samples, sample)
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
  samples = c()
  
  for (file in unpair(allfiles)) {
    
    if (paired==T) { sample = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file)) 
    } else { sample = gsub(".fastq.gz$","",basename(file)) }
    samples = c(samples, sample)
    num = as.numeric( system(sprintf("zcat %s | wc -l | awk '{print $1/4}'", file), intern=T) )
    num_reads = c( num_reads, num )  
    cat(sample, num, "\n")
  }
  report = data.frame(samples=samples, FASTQ=num_reads)
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
  samples = c()
  counts = c()
  for (file in unpair(files)) {
    if (paired==T) { sample = gsub("_R1.fastq.gz$|_R2.fastq.gz$","",basename(file)) 
    } else { sample = gsub(".fastq.gz$","",basename(file)) }
    samples = c( samples, sample )
    count = as.numeric( system(sprintf("zcat %s | wc -l | awk '{print $1/4}'", file), intern=T) )
    counts = c( counts, count )
    cat(paste(sample, count, "\n"), sep=" ")
  }
  df_counts = data.frame(samples=samples, TrimGalore=counts)
  report = merge(report, df_counts, by="samples")
  report$TrimGalore_perc = report$TrimGalore / report$FASTQ * 100
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
  samples = c()
  
  for (file in files) {
    sample = gsub(".bam$","",basename(file))
    samples = c(samples, sample)

    ## count mapped (-F4) and primary (-F256) fragments
    count = as.numeric( system(sprintf("%s view -c -F260 %s ", file.path(samtools_dir,"samtools"), file), intern=T) )
    if (paired==T) { count = round(count / 2) }

    counts = c( counts, count )
    cat(paste(sample, count, "\n"), sep=" ")
  }
  df_counts = data.frame(samples=samples, HISAT2=counts)
  report = merge(report, df_counts, by="samples")
  report$HISAT2_perc = report$HISAT2 / report$FASTQ * 100
  report
}


## TopHat2
indir = file.path(main_outdir,"TopHat2")
if ( file.exists( indir ) ){
  files <- sort(list.files(indir, pattern="*.bam$", full.names=T))
  ## files = sample(files) ## debugging only!!!
  files
  counts = c()
  samples = c()
  
  for (file in files) {
    sample = gsub(".bam$","",basename(file))
    samples = c(samples, sample)
    
    ## count mapped (-F4) and primary (-F256) fragments
    count = as.numeric( system(sprintf("%s view -c -F260 %s ", file.path(samtools_dir,"samtools"), file), intern=T) )
    if (paired==T) { count = round(count / 2) }

    counts = c( counts, count )
    cat(paste(sample, count, "\n"), sep=" ")
  }
  df_counts = data.frame(samples=samples, TopHat2=counts)
  report = merge(report, df_counts, by="samples")
  report$TopHat2_perc = report$TopHat2 / report$FASTQ * 100
  report
}


## STAR
indir = file.path(main_outdir,"STAR")
if ( file.exists( indir ) ){
  files <- sort(list.files(indir, pattern="*.bam$", full.names=T))
  ## files = sample(files) ## debugging only!!!
  files
  counts = c()
  samples = c()

  for (file in files) {
    sample = gsub(".bam$","",basename(file))
    samples = c(samples, sample)

    ## count mapped (-F4) and primary (-F256) fragments
    count = as.numeric( system(sprintf("%s view -c -F260 %s ", file.path(samtools_dir,"samtools"), file), intern=T) )
    if (paired==T) { count = round(count / 2) }

    counts = c( counts, count )
    cat(paste(sample, count, "\n"), sep=" ")
  }
  df_counts = data.frame(samples=samples, STAR=counts)
  report = merge(report, df_counts, by="samples")
  report$STAR_perc = report$STAR / report$FASTQ * 100
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

 df_counts = data.frame( row.names=NULL,samples=colnames(counts), featureCounts=apply(counts,2,sum) )
 report = merge(report, df_counts, by="samples")
 report$featureCounts_perc = report$featureCounts / report$FASTQ * 100
 report
}

##htseq-count
indir = file.path(main_outdir,"htseq-count")
infile = file.path(indir,"counts.txt")
if ( file.exists( infile ) ) {
  cat(infile, "\n") 
  counts = read.table(infile, header=TRUE)
  head(counts)
  
  df_counts = data.frame( row.names=NULL, samples=colnames(counts), htseq_count=apply(counts,2,sum) )
  report = merge(report, df_counts, by="samples")
  report$htseq_count_perc = report$htseq_count / report$FASTQ * 100
  report
}

rownames(report) = report$samples
report$samples = NULL
write.table(report,"Report.tsv", sep="\t", quote=FALSE, col.names=NA) # save to file
