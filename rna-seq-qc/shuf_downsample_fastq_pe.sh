#!/usr/bin/bash

## Usage: bash shuf_downsample_fastq_se.sh 1000 a.fastq.gz b.fastq.gz a.downsampled.fastq.gz b.downsampled.fastq.gz
## Downsampling has NO seed!!!

DIFF=$(diff <(zcat ${2} | sed -n '1p;1~4p' | cut -d " " -f1) <(zcat ${3} | sed -n '1p;1~4p' | cut -d " " -f1)) 
if [ "$DIFF" != "" ] 
then
    echo "Error! FASTQ sequence identifiers do NOT match!"
    exit 1
else
    paste <(zcat ${2} | sed '/^$/d' | paste -d "\t" - - - - ) <(zcat ${3} | sed '/^$/d' | paste -d "\t" - - - - ) \
    | shuf \
    | head -${1} \
    | tee >(cut -f1,2,3,4 | sed 's/\t/\n/g' | gzip > ${4}) >(cut -f5,6,7,8 | sed 's/\t/\n/g' | gzip > ${5}) >/dev/null
fi
