#!/usr/bin/env bash

## reservoir sampling
awk -v k=$1 '{ s=x++<k?x-1:int(rand()*x); if (s<k) R[s]=$0 } END { for (i in R) print R[i] }' $2
