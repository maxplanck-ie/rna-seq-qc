#!/usr/bin/env python

import sys
import random

if len(sys.argv) == 3:
    input = open(sys.argv[2],'r')
elif len(sys.argv) == 2:
    input = sys.stdin;
else:
    sys.exit("Usage:  python reservoir.py <lines> <?file>")

N = int(sys.argv[1]);
sample = [];

for i,line in enumerate(input):
    if i < N:
        sample.append(line)
    elif i >= N and random.random() < N/float(i+1):
        sample[random.randint(0,len(sample)-1)] = line

for line in sample:
    sys.stdout.write(line)
