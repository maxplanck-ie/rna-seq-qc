library(gdata)

.Library
.Library.site
.libPaths()

sessionInfo()

args = commandArgs(TRUE)

## For debugging only!!! #######################################################
# setwd("/data/processing/kilpert/test/rna-seq-qc/Ausma/PE_mm10_FULL_hisat_trim/project_report/")
# args = c('/data/manke/kilpert/datasets/Ausma/',
#         '/data/processing/kilpert/test/rna-seq-qc/Ausma/PE_mm10_FULL_hisat_trim/')
################################################################################

main_indir = args[1]
main_outdir = args[2]


## binaries ####################################################################
hostname = Sys.info()["nodename"]
if (hostname == "pc305.immunbio.mpg.de") {
  samtools_dir = "/home/kilpert/Software/samtools/samtools-1.2/"
} else if ( startsWith(hostname,'deep') ) {
  samtools_dir = "/package/samtools-1.2/bin/"
} else {
  samtools_dir = ""
}


print("Generating project report...")


## get R1 read if paired ####################################################

unpair <- function(allfiles){
  files = c()
  names = c()
  for (file in allfiles) {
    name = gsub("_R1+.fastq.gz$|_R2+.fastq.gz$","",basename(file))
    
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
    name = gsub("_R1+.fastq.gz$|_R2+.fastq.gz$","",basename(file))
    
    if (!name %in% names) {
      names = c(names, name)
    } else {
      paired = T
    }
  }  
  return(paired)
}



## FASTQ: total ################################################################
indir = main_indir
if (file.exists(indir)){
  allfiles <- list.files(indir, pattern="*.fastq.gz$", full.names=T)
  allfiles
  paired = is_paired(allfiles)
  
  num_reads = c()
  files = c()
  names = c()
  
  for (file in unpair(allfiles)) {
    name = gsub("_R1+.fastq.gz$|_R2+.fastq.gz$","",basename(file))
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
  files <- list.files(indir, pattern="*.fastq.gz$", full.names=T)
  files
  names = c()
  num_reads = c()
  for (file in unpair(files)) {
    name = gsub("_R1+.fastq.gz$|_R2+.fastq.gz$","",basename(file))
    names = c( names, name )
    num = as.numeric( system(sprintf("zcat %s | wc -l | awk '{print $1/4}'", file), intern=T) )
    num_reads = c( num_reads, num )
    cat(paste(name, num, "\n"), sep=" ")
  }
  report$TrimGalore = num_reads
  report$TrimGalore_perc = apply(report, 1, function(x) x["TrimGalore"]/x["FASTQ"]*100)
  report
} 


## BAM: mapped #################################################################
## HISAT
indir = file.path(main_outdir,"HISAT")
if ( file.exists( indir ) ){
  files <- list.files(indir, pattern="*.bam$", full.names=T)
  files
  num_reads = c()
  
  for (file in files) {
    name = gsub(".bam$","",basename(file))
    
    if (paired==T) {
      num = as.numeric( system(sprintf("%s view -c -f3 -F256 %s | awk '{print $1/2}'", file.path(samtools_dir,"samtools"), file), intern=T) )
    } else {
      num = as.numeric( system(sprintf("%s view -c -F256 %s | awk '{print $1/2}'", file.path(samtools_dir,"samtools"), file), intern=T) )
    }
    num_reads = c( num_reads, num )
    cat(paste(name, num, "\n"), sep=" ")
  }
  report$HISAT = num_reads
  report$HISAT_perc = apply(report, 1, function(x) x["HISAT"]/x["FASTQ"]*100)
  report
}

## TopHat2
indir = file.path(main_outdir,"TopHat2")
if ( file.exists( indir ) ){
  files <- list.files(indir, pattern="*.bam$", full.names=T)
  files
  num_reads = c()
  
  for (file in files) {
    name = gsub(".bam$","",basename(file))
    
    if (paired==T) {
      num = as.numeric( system(sprintf("%s view -c -f3 -F256 %s | awk '{print $1/2}'", file.path(samtools_dir,"samtools"), file), intern=T) )
    } else {
      num = as.numeric( system(sprintf("%s view -c -F256 %s | awk '{print $1/2}'", file.path(samtools_dir,"samtools"), file), intern=T) )
    }
    num_reads = c( num_reads, num )
    cat(paste(name, num, "\n"), sep=" ")
  }
  report$TopHat2 = num_reads
  report$TopHat2_perc  = apply(report, 1, function(x) x["TopHat2"]/x["FASTQ"]*100)
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
 report$featureCounts = apply(counts, 2, sum)
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
  report$htseq-count = apply(counts, 2, sum)
  report$htseq-counts_perc = apply(report, 1, function(x) x["htseq-count"]/x["FASTQ"]*100)
  report
}

write.table(report,"Report.tsv", sep="\t", quote=FALSE, col.names=NA) # save to file

