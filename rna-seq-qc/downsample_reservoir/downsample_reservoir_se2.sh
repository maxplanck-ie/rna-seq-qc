#!/usr/bin/bash

## Usage: bash reservoir_downsample_se.sh 1000 a.fastq.gz a.downsampled.fastq.gz
## Downsampling has NO seed!!!

bindir=$(dirname $0)

pigz -dc ${2} | sed '/^$/d' | paste -d "\t" - - - - | $bindir/reservoir ${1} | sed 's/\t/\n/g' | pigz -9 > ${3}
