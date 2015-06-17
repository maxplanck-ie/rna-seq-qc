rna-seq-qc v0.6.1
=================

RNA-seq pipeline for processing RNA sequence data from high throughput sequencing.

Fabian Kilpert - June 17, 2015
email: kilpert@ie-freiburg.mpg.de

This software is distributed WITHOUT ANY WARRANTY!
--------------------------------------------------

Following steps are executed in succession: FASTQ subsampling (optiona), quality check
with FASTQC, trimming of reads with Trim Galore (optional), estimation of insert size
and strand specificity with RSeQC, mapping with TopHat, extensive quality check with RSeQC,
counting of features with featureCounts (default) or htseq-count, differential expression
analysis with DESeq2 (optional).

The pipeline requires gzipped FASTQ files (.fastq.gz) for input, which are loaded
from an input directory (-i INDIR). Read files belonging together require the exact same
base name but ending either in "_R1" or "_R2" right in front of the .fastq.gz extension.
(e.g. reads_R1.fastq.gz, reads_R2.fastq.gz). In addition, a specific genome version argument
must be provided (e.g. -g mm10) to define the reference data used for annotation.
This loads a number of indexes for mapping programs (Bowtie2, TopHat2, HISAT etc.) from the
corresponding configuration file of the rna-seq-qc sub-folder (e.g. rna-seq-qc/mm10.cfg).
Additional genomes for selection can be provided as cfg-files by the user. The pipeline
works for single end and paired end sequences alike.

The DE analysis is only executed if a valid setup table is provided
(e.g. --DE setup_table.tsv), which defines the relationships of the samples.

More information on the pipeline can be found on the wiki page:
http://epicenter/wiki/index.php/RNA-seq_pipeline_(rna-seq-pc.py)

Example:
    python rna-seq-qc.py -i /path/to/fastq_dir -o /path/to/ouput_dir -g mm10 -v --DE sampleInfo.tsv