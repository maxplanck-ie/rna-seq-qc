#!/usr/local/bin/python

__description__ = """
    downsample_fastq v0.2
    Fabian Kilpert - February 9, 2014
    email: kilpert@ie-freiburg.mpg.de
    ---------------------------------
    Down-sampling, random or non-random, of a FASTQ file.

    This software is distributed WITHOUT ANY WARRANTY!

    Example:
        python downsample_fastq.py -v -n 200000 -s 12345 example.fastq.gz few.fastq.gz
    """


import argparse
import datetime
import gzip
import random
import shutil
import subprocess
import sys
import textwrap

# import socket
# if socket.gethostname() == "pc196":
#     sys.argv = [sys.argv[0],
#                 '-v',
#                 '-s','123',
#                 '-n','1000',
#                 #'/data/processing/heyne/data_miRNA_knockdown_fabian/KDcl2_R1.fastq.gz',
#                 '/data/processing/kilpert/downsampled/KDcl2_R1.10.fastq.gz',
#                 '/data/processing/kilpert/downsampled/KDcl2_R1.fastq.gz'
#                 ]

#### COMMAND LINE PARSING ##############################################

def parse_args():
    """
    Parse arguments from the command line.
    """
    parser = argparse.ArgumentParser(
                        prog='downsample_fastq.py',
                        formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(__description__))

    parser.add_argument("-n", "--number", dest="num", metavar="INT", help="Number of reads (default: 100)", type=int, default=100)
    parser.add_argument("-s", "--seed", dest="seed", metavar="INT", help="Random number seed", type=int, default=None)
    parser.add_argument("--head", dest="head", action="store_true", default=False, help="Down-sampling (non-random!) just from the head of the file (default: False)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose output")

    ### Positional arguments (required!)
    parser.add_argument("infile", type=str, help="Input file")
    parser.add_argument("outfile", type=str, help="Output file")

    args = parser.parse_args()
    
    args.rnd = not args.head
    return args


#### COUNTING LINES IN FILE ############################################

def count_lines(infile):
    """
    Count the reads in a FASTQ file.
    """
    return int(subprocess.check_output("zcat {} | wc -l".format(infile), shell=True))


def count_lines2(infile):
    """
    Count the reads in a FASTQ file.
    """
    return int(sum(1 for line in gzip.open(infile, "r")))


def count_lines3(infile):
    f = gzip.open(infile)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines


#### DOWNSAMPLING ######################################################

def downsample_head(infile, outfile, n=100, verbose=False):
    """
    Downsample n reads from the head of the file.
    """
    if verbose:
        print "Non-random downsampling (from head of the file)"
    with gzip.open(infile, "r") as f1:
        with gzip.open(outfile, "w") as f2:
            f2.writelines(f1.readline() for i in xrange(n*4))


def downsample_random(infile, outfile, n=100, seed=None, verbose=False):
    """
    Downsample randomly from FASTQ file.
    """
    if verbose:
        print "Random downsampling"
        print "Seed:", seed

    total = count_lines3(infile)/4
    if verbose:
        print "Total:", "{:,}".format(total)

    if total <= n:
        shutil.copyfile(infile, outfile)
        exit(0)

    if seed:
        random.seed(seed)
    else:
        random.seed()
    read_list = sorted(random.sample(xrange(total), n), reverse=True)
    if verbose:
        print "Requested:", "{:,}".format(len(read_list))

    with gzip.open(infile, "r") as f1:
        with gzip.open(outfile, "w") as f2:
            for i in range(total):
                four_lines = [f1.readline() for x in range(4)]

                if i and i % 1000000 == 0:
                    if verbose:
                        print "{:,} : {:,} ({:.0%})".format(i, n-len(read_list), (1.0*n-len(read_list))/n )
                    else:
                        sys.stdout.write( "\r{:.0%}".format( (1.0*n-len(read_list))/n ) )
                        sys.stdout.flush()

                if i == read_list[-1]:
                    f2.writelines(four_lines)
                    read_list.pop()
                    if not read_list:
                        if verbose:
                            print "{:,} : {:,} ({:.0%})".format(i, n-len(read_list), (1.0*n-len(read_list))/n )
                        else:
                            sys.stdout.write( "\r{:.0%}".format( (1.0*n-len(read_list))/n ) )
                            sys.stdout.flush()
                        break


#### MAIN ##############################################################

def main():
    """
    Main program.
    """
    args = parse_args()
    #print "Args:", args

    t1 = datetime.datetime.now()
    
    if not args.rnd:
        downsample_head(args.infile, args.outfile, n=args.num, verbose=args.verbose)
    else:
        downsample_random(args.infile, args.outfile, n=args.num, seed=args.seed, verbose=args.verbose)

    t2 = datetime.datetime.now()
    if args.verbose:
        print "Duration:", t2 - t1


if __name__ == "__main__":
    main()
