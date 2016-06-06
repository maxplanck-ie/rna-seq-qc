#!/usr/bin/bash

## Usage: bash reservoir_downsample_se.sh 1000 a.fastq.gz a.downsampled.fastq.gz
## Downsampling has NO seed!!!

bindir=$(dirname $0)

zcat ${2} | sed '/^$/d' | sed 's/\t/ /' | paste -d "\t" - - - - | $bindir/reservoir ${1} | sed 's/\t/\n/g' | gzip > ${3}
