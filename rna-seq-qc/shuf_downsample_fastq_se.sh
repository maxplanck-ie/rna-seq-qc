#!/usr/bin/bash

## Usage: bash shuf_downsample_fastq_se.sh 1000 a.fastq.gz a.downsampled.fastq.gz
## Downsampling has NO seed!!!

zcat ${2} | sed '/^$/d' | paste -d "\t" - - - - | shuf | head -${1} | sed 's/\t/\n/g' | gzip > ${3}
