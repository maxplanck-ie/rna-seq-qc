#!/usr/local/bin/python

__description__ = """
    rna-seq-qc v0.12
    Fabian Kilpert, 10-13-2014
    email: kilpert@ie-freiburg.mpg.de
    ---------------------------------
    RNA-seq pipeline for processing RNA sequence data from high throughput sequencing.

    This software is distributed WITHOUT ANY WARRANTY!

    Following steps are executed in succession: FASTQ subsampling (optiona), quality check
    with FASTQC, trimming of reads with Trim Galore (optional), estimation of insert size
    and strand specificity with RSeQC, mapping with TopHat, extensive quality check with RSeQC,
    counting of features with featureCounts (default) or htseq-count, differential expression
    analysis with DESeq2 (optional).

    The pipeline requires gzipped FASTQ files (.fastq.gz) for input, which are loaded
    from an input directory (-d INDIR). Read files belonging together require the exact same
    base name but ending either in "_R1" or "_R2" right in front of the .fastq.gz extension.
    (e.g. reads_R1.fastq.gz, reads_R2.fastq.gz). In addition, a specific genome version argument
    must be provided (e.g. -g mm10) to define the reference data used for annotation.
    This loads a number of indexes for mapping programs (Bowtie2, TopHat2, etc.) from the
    corresponding configuration file of the rna-seq-qc sub-folder (e.g. rna-seq-qc/mm10.cfg).
    Additional genomes for selection can be provided as cfg-files by the user. The pipeline
    works for single end and paired end sequences alike.

    The DE analysis is only executed if a valid setup table is provided
    (e.g. --DE setup_table.tsv), which defines the relationships of the samples.

    More information on the pipeline can be found on the wiki page:
    http://epicenter/wiki/index.php/RNA-seq_pipeline_(rna-seq-pc.py)

    Example:
        python rna-seq-qc.py -d /path/to/fastq_dir -o /path/to/ouput_dir -g mm10 -v --DE sampleInfo.tsv
    """

import argparse
from collections import OrderedDict
import datetime
import gzip
import os
import os.path
from Queue import Queue
import random
import re
import shutil
import socket
import subprocess
import sys
import textwrap
from threading import Thread
import time

########################################################################################################################
## for pc196 (developer) only!
########################################################################################################################
#print socket.gethostname()
if socket.gethostname() == "pc196":

    # 2_samples
    # sys.argv = [sys.argv[0],
    #             '-d', '/data/projects_2/kilpert/datasets/Bottomly/2_samples/',
    #             '-o', '/data/projects_2/kilpert/dev/rna-seq-pipeline/2_samples',
    #             '-g','mm10',
    #             '-v',]

    # ## 140731_MiSeq_Ausma
    # sys.argv = [sys.argv[0],
    #             '-d', '/data/projects_2/kilpert/dev/rna-seq-pipeline/test_data/140731_MiSeq_Ausma',
    #             '-o', '140731_MiSeq_Ausma',
    #             '-g', 'mm10',
    #             '-v',
    #             '--DE', '/data/projects_2/kilpert/dev/rna-seq-pipeline/sampleInfo.tsv'
    #             ]

    ## SE500
    # sys.argv = [sys.argv[0],
    #             '-d', '/data/projects_2/kilpert/dev/rna-seq-pipeline/test_data/SE',
    #             '-o', '/data/projects_2/kilpert/dev/rna-seq-pipeline/outdir500_SE',
    #             '--fastq-downsample', '400',
    #             '-g', 'mm10',
    #             #'--overwrite',
    #             '-v',
    #             #'--no-bam',
    #             #'--no-trim',
    #             #'--trim_galore', '-q 30'
    #             '--count-prg', 'featureCounts'
    #             ]

    # ## PE500
    # sys.argv = [sys.argv[0],
    #             '-d', '/data/projects_2/kilpert/dev/rna-seq-pipeline/test_data/PE',
    #             '-o', '/data/projects_2/kilpert/dev/rna-seq-pipeline/outdir500_PE',
    #             '--fastq-downsample', '400',
    #             '-g', 'mm10',
    #             #'--overwrite',
    #             '-v',
    #             #'--no-bam',
    #             #'--no-trim',
    #             '--trim_galore', '-q 30',
    #             #'--insert-metrics','Picard'
    #             '--count-prg', 'featureCounts'
    #             ]

    ## PE full
    # #/home/kilpert/git/rna-seq-pipeline/rna-seq-qc10.py -d /data/projects_2/kilpert/megumi/02_new/00_data -o /data/projects_2/kilpert/dev/rna-seq-pipeline/outdir_PE_full -g mm10 -v
    # sys.argv = [sys.argv[0],
    #             '-d', '/data/projects_2/kilpert/megumi/02_new/00_data',
    #             '-o', '/data/projects_2/kilpert/dev/rna-seq-pipeline/outdir_PE_full',
    #             '-g', 'mm10',
    #             '-v'
    #             ]

    # ## SE full
    # sys.argv = [sys.argv[0],
    #             '-d', '/data/projects/kilpert/datasets/Bottomly',
    #             '-o', '/data/projects_2/kilpert/dev/rna-seq-pipeline/outdir_SE_full',
    #             '-g', 'mm10',
    #             '-v',
    #             ]

    if "--trim_galore" in sys.argv:
        sys.argv2 = sys.argv
        sys.argv2[sys.argv.index( '--trim_galore')+1] = "'" + sys.argv[sys.argv2.index( '--trim_galore')+1] + "'"   # just to output all trim galore option as one string
        print " ".join(sys.argv2)   # output command line
    if "--featureCounts" in sys.argv:
        sys.argv2 = sys.argv
        sys.argv2[sys.argv.index( '--featureCounts')+1] = "'" + sys.argv[sys.argv2.index( '--featureCounts')+1] + "'"   # just to output all trim galore option as one string
        print " ".join(sys.argv2)   # output command line
    else:
        print " ".join(sys.argv)   # output command line


#### PATHS #############################################################################################################
## Path defaults to all needed scripts and programs

bowtie2_path = "bowtie2"            #version 2.2.3
fastqc_path = "/package/FastQC_v0.11.2/fastqc"
htseq_count_path = "/package/HTSeq-0.6.1/bin/htseq-count"
picardtools_dir_path = "/package/picard-tools-1.121/"
rseqc_dir_path = "/package/RSeQC-2.4/bin/"
samtools_path = "/package/samtools-1.1/samtools"
tophat2_path = "/package/tophat-2.0.13.Linux_x86_64/tophat2"
trim_galore_path = "/package/trim_galore_v0.3.7/trim_galore"
ucsctools_dir_path = "/package/UCSCtools/"
feature_counts_path = "/package/subread-1.4.5-p1-Linux-x86_64/bin/featureCounts"
deseq2_path = "rna-seq-qc/DESeq2.R"     # relative to main script

## Different configurations
if socket.gethostname() == "pc196":
    fastqc_path = "fastqc"
    htseq_count_path = "htseq-count"
    samtools_path = "samtools"
    tophat2_path = "tophat"
    rseqc_dir_path = ""
    ucsctools_dir_path = ""


#### DEFAULT VARIABLES #################################################################################################
default_processors = 2          # Total number of processors shared by all threads
default_threads = 1             # Number of threads for parallel file processing
samtools_max_mem = 4
if socket.gethostname() == "deep1":
    default_processors = 20
    default_threads = 2
    samtools_max_mem = 20
is_error = False

## Software variables
rseqc2tophat = {"1++,1--,2+-,2-+": "fr-secondstrand", "1+-,1-+,2++,2--": "fr-firststrand",
                "++,--": "fr-secondstrand", "+-,-+": "fr-firststrand"}
rseqc2htseq = {"1++,1--,2+-,2-+": "yes", "1+-,1-+,2++,2--": "reverse",
                "++,--": "yes", "+-,-+": "reverse"}

downsample_size = 1000000       #for CollectInsertSizeMetrics


#### COMMAND LINE PARSING ##############################################################################################

def parse_args():
    """
    Parse arguments from the command line.
    """
    parser = argparse.ArgumentParser(
    prog='rna-seq-qc.py',
    formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(__description__))

    ### Optional arguments
    parser.add_argument('--version', action='version', version='%(prog)s v0.12')
    parser.add_argument("-d", "--input-dir", dest="indir", required=True, help="Input dir containing (FASTQ)")
    parser.add_argument("-o", "--output-dir", dest="outdir", required=True, help="Output directory")
    parser.add_argument("--overwrite", dest="overwrite", action = "store_true", default=False, help="Overwrite results in existing folders!")
    parser.add_argument("-p", "--processors", dest="processors", metavar="INT", help="Maximum number of processors available shared by all threads (default: {}).format(default_processors)", type=int, default=default_processors)
    parser.add_argument("--threads", dest="threads", metavar="INT", help="Number of threads for running files in parallel (default: {})".format(default_threads), type=int, default=default_threads)
    parser.add_argument("--seed", dest="seed", metavar="INT", help="Random number seed", type=int, default=None)
    parser.add_argument("--fastq-downsample", dest="fastq_downsample", metavar="INT", help="Subsample first n fastq sequences", type=int, default=None)
    parser.add_argument("-g", "--genome", dest="genome", required=True, help="Reference genome build")
    parser.add_argument("--trim_galore", dest="trim_galore_opts", metavar="STR", help="Trim Galore! option string (default: '--stringency 2')", type=str, default="--stringency 2")
    parser.add_argument("--tophat_opts", dest="tophat_opts", metavar="STR", help="TopHat2 option string", type=str, default="")     #--library-type fr-firststrand
    parser.add_argument("--bowtie_opts", dest="bowtie_opts", metavar="STR", help="Bowtie2 option string (default: '--end-to-end --fast')", type=str, default="--end-to-end --fast")
    parser.add_argument("--featureCounts_opts", dest="featureCounts_opts", metavar="STR", help="featureCounts option string (default: '')", type=str, default="-Q 10")
    parser.add_argument("--htseq-count_opts", dest="htseq_count_opts", metavar="STR", help="HTSeq htseq-count option string", type=str, default="--mode union")
    parser.add_argument("--insert-metrics", dest="insert_metrics", metavar="STR", help="Calculate insert size metrics (mean, sd) using RSeQC (default) or Picard", type=str, default="RSeQC")
    parser.add_argument("--count-prg", dest="count_prg", metavar="STR", help="Program used for counting features: htseq-count (default) or featureCounts", type=str, default="featureCounts")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose output")
    parser.add_argument("--no-bam", dest="no_bam", action="store_true", default=False, help="First steps only. No alignment. No BAM file.")
    parser.add_argument("--no-trim", dest="no_trim", action="store_true", default=False, help="Do not trim FASTQ reads. Default: Use Trim Galore! with default settings.")
    parser.add_argument("--DE", dest="sample_info", help="Information on samples (required for DE analysis)")

    ### Positional arguments (required!)
    args = parser.parse_args()

    ### Add useful paths
    args.script = os.path.dirname(os.path.realpath(__file__))
    args.cwd = os.getcwd()
    args.now = datetime.datetime.now()

    ### Sanity checking
    ## Input dir
    ##print "indir:", args.indir
    args.indir = os.path.expanduser(args.indir)
    args.indir = os.path.realpath(args.indir)
    if (not os.path.isdir(args.indir)):
        print "Error: The specified input directory does not exist or is not a directory:\n", args.indir
        exit(1)

    if not args.sample_info is None:
        args.sample_info = os.path.expanduser(args.sample_info)
        args.sample_info = os.path.realpath(args.sample_info)
        if (not os.path.isfile(args.sample_info)):
            print "Error: The specified file does not exist or is not readable:\n", args.sample_info
            exit(1)

    ## Output dir
    args.outdir = os.path.join(args.cwd, os.path.expanduser(args.outdir))

    args.error = 0

    ### Add
    if args.seed is None:
        args.seed = random.randint(12345, 987654321)

    ## Get reference data paths from config file
    ref_cfg_file_path = os.path.join(args.script, "rna-seq-qc/{}.cfg".format(args.genome))
    if not os.path.isfile(ref_cfg_file_path):
        print "Error! Configuration file NOT found for {}: {}".format(args.genome, ref_cfg_file_path)
        exit(1)
    configs = parse_config(ref_cfg_file_path)
    try:
        args.fasta_index = configs["fasta_index"]
        args.bowtie2_index = configs["bowtie2_index"]
        args.tophat2_index = configs["tophat2_index"]
        args.gtf = configs["gtf"]
        args.bed = configs["bed"]
    except:
        print "Error! Unable to read paths from config file:", ref_cfg_file_path
        exit(1)

    ## BioMart
    args.biomart = os.path.join(args.script, "rna-seq-qc/BioMart.{}.tsv".format(args.genome))
    if not os.path.isfile(args.biomart):
        args.biomart = ""
        
    return args


#### GENERIC FUNCTIONS #################################################################################################

class Qjob():
    def __init__(self, cmds=None, logfile=None, shell=False):
        self.cmds = cmds
        self.logfile = logfile
        self.shell = shell
    def run(self):
        pass


def parse_config(file_path):
    """
    Parse a configuration text file, e.g. *.cfg
    """
    options = {}
    for line in file(file_path):
        line = line.rstrip("\r\n")
        if not line.startswith("#"):
            if "=" in line:
                key, value = line.split('=', 1)
                value = value.strip()
                value = value.strip('"')
                value = value.strip("'")
                options[key] = value
    return options


def check_for_paired_infiles(args, indir, ext, verbose=False):
    """
    Read input files from given dir that have the user specified extension.
    """
    infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(ext)])

    paired_infiles = OrderedDict()
    for infile in infiles:
        fname = os.path.basename(infile).replace(ext, "")
        m = re.match("^(\w+)_R*[1|2]$", fname)
        if m:
            bname = m.group(1)
            if bname not in paired_infiles:
                paired_infiles[bname] = [infile]
            else:
                paired_infiles[bname].append(infile)
    if paired_infiles:
        for pair in paired_infiles.values():
            if len(pair) != 2:
                print "Error! Unpaired input file:", pair
                exit(1)
        infiles = paired_infiles.values()

    ## Paired
    if len(infiles[0]) == 2:
        args.paired = True
        print "Paired end FASTQ files ({} pairs)".format(len(infiles))
        if verbose:
            for pair in infiles:
                print "  {}  {}".format(os.path.basename(pair[0]), os.path.basename(pair[1]))
    else:
    ## Single
        args.paired = False
        print "Single end FASTQ files ({})".format(len(infiles))
        if verbose:
            for infile in infiles:
                print "  {}".format(os.path.basename(infile))

    return infiles


def add_leading_zeros(num, length=2):
    """
    Add leading zeros to number.
    """
    padding = length - len("{}".format(num))
    return "{}{}".format(padding*"0", num)


def run_subprocess(cmd, shell=False, logfile=None, verbose=False):
    """
    Run the subprocess.
    """
    if shell == True:
        return subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    else:
        if type(cmd) is str:
            cmd = cmd.split()
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if verbose:
                print "PID:", p.pid
            (stdoutdata, stderrdata) = p.communicate()
            if logfile:
                with open(logfile, "a+") as log:
                    if stderrdata:
                        log.write("{}\n".format(stderrdata))
                    if stdoutdata:
                        log.write("{}\n".format(stdoutdata))
            if stderrdata:
                print stderrdata
            if stdoutdata:
                print stdoutdata
            if p.returncode != 0:
                print "Error! Command failed with return code:", p.returncode
                print " ".join(cmd) + "\n"
        finally:
            if p.poll() == None:
                print "Error! Trying to terminate PID {}".format(p.pid)
                p.terminate()
                time.sleep(20)
                if p.poll() == None:
                    print "Error! Trying to kill PID {}!".format(p.pid)
                    p.kill()
                exit(1)
        return p.wait()


def queue_worker(q, verbose=False):
    """
    Worker executing jobs (command lines) sent to the queue.
    """
    global is_error

    while True:
        job = q.get()
        if job.logfile:
            logfile = job.logfile
        else:
            logfile = None
        for cmd in job.cmds:
            if job.logfile:
                with open(job.logfile, "a+") as log:
                    log.write("{}\n\n".format(cmd))
            if verbose:
                print cmd
            return_code = run_subprocess(cmd, job.shell, logfile, verbose)
            if return_code:
                is_error = True
        q.task_done()


def get_strand_from_rseqc(infile, prog):
    """
    Infer strand-specificity.
    """
    specificity_list = []
    is_first_file = True

    with open(infile, "r") as f:
        strands = OrderedDict()
        for line in f.readlines():
            line = line.strip()
            if line.startswith('Fraction of reads explained by "'):
                values = line.replace('Fraction of reads explained by "','').split('": ')
                strands[values[0]] = float(values[1])
            if line.startswith('Fraction of reads explained by other combinations: '):
                value = float(line.replace('Fraction of reads explained by other combinations: ',''))
                if value >= 0.2:
                    print "Error! A larger fraction ({}) of reads is explained by uncommon strand combinations!".format(value)
                    print infile
                    exit(1)

        if len(strands.keys()) != 2:
            print "Error! Unclear strand-specificity in:"
            print infile
            exit(1)

        threshold = 0.8      #min quotient threshold

        k = strands.keys()
        v = strands.values()

        if prog == "rseqc":
            specificity = None
            if '++,--' in strands.keys() and '+-,-+' in strands.keys():
                if strands['++,--'] >= threshold and strands['+-,-+'] <= threshold:
                    specificity = '++,--'
                elif strands['++,--'] <= threshold and strands['+-,-+'] >= threshold:
                    specificity = '+-,-+'
            if '1++,1--,2+-,2-+' in strands.keys() and '1+-,1-+,2++,2--' in strands.keys():
                if strands['1++,1--,2+-,2-+'] >= threshold and strands['1+-,1-+,2++,2--'] <= threshold:
                    specificity = '1++,1--,2+-,2-+'
                elif strands['1++,1--,2+-,2-+'] <= threshold and strands['1+-,1-+,2++,2--'] >= threshold:
                    specificity = '1+-,1-+,2++,2--'

        elif prog == "tophat":
            specificity = "fr-unstranded"
            if '++,--' in strands.keys() and '+-,-+' in strands.keys():
                if strands['++,--'] >= threshold and strands['+-,-+'] <= threshold:
                    specificity = rseqc2tophat['++,--']             #fr-secondstrand
                elif strands['++,--'] <= threshold and strands['+-,-+'] >= threshold:
                    specificity = rseqc2tophat['+-,-+']             #fr-firststrand
            if '1++,1--,2+-,2-+' in strands.keys() and '1+-,1-+,2++,2--' in strands.keys():
                if strands['1++,1--,2+-,2-+'] >= threshold and strands['1+-,1-+,2++,2--'] <= threshold:
                    specificity = rseqc2tophat['1++,1--,2+-,2-+']   #fr-secondstrand
                elif strands['1++,1--,2+-,2-+'] <= threshold and strands['1+-,1-+,2++,2--'] >= threshold:
                    specificity = rseqc2tophat['1+-,1-+,2++,2--']   #fr-firststrand

        elif prog == "htseq":
            specificity = "no"
            if '++,--' in strands.keys() and '+-,-+' in strands.keys():
                if strands['++,--'] >= threshold and strands['+-,-+'] <= threshold:
                    specificity = rseqc2htseq['++,--']              #yes
                elif strands['++,--'] <= threshold and strands['+-,-+'] >= threshold:
                    specificity = rseqc2htseq['+-,-+']              #reverse
            if '1++,1--,2+-,2-+' in strands.keys() and '1+-,1-+,2++,2--' in strands.keys():
                if strands['1++,1--,2+-,2-+'] >= threshold and strands['1+-,1-+,2++,2--'] <= threshold:
                    specificity = rseqc2htseq['1++,1--,2+-,2-+']    #yes
                elif strands['1++,1--,2+-,2-+'] <= threshold and strands['1+-,1-+,2++,2--'] >= threshold:
                    specificity = rseqc2htseq['1+-,1-+,2++,2--']    #reverse

    if is_first_file:
        specificity_list.append(specificity)
    else:
        is_first_file = False
        if specificity not in specificity_list:
            print "Error! Multiple strand specificities detected within different samples."
            print " ".join(specificity_list)
            exit(1)
        else:
            specificity_list.append(specificity)

    return specificity_list[0]


def downsample_fastq(infile, outfile, number=1000000, rnd=False, seed=None):
    """
    Downsample FASTQ file. Non-random is fast and samples just the from top of the input file.
    RND mode reads the full input file into memory and applies a pseudo random selection.
    """
    if not rnd:
        with gzip.open(infile, "r") as f1:
            with gzip.open(os.path.basename(outfile), "w") as f2:
                f2.writelines(f1.readline() for i in xrange(number*4))      # Just grab the first n*4 lines from input files!!!
    else:
        lines = []
        with gzip.open(infile, "r") as f:
            for line in f:
                line = line.strip()
                lines.append(line)
        ##print "FASTQ reads:", len(lines)/4
        fastq = []
        for i in xrange(0, len(lines),4):
            if not lines[i].startswith("@"):
                print >> sys.stderr, "Error: Unusual SEQ-ID in line:" + i + "\n", lines[i]
                exit(1)
            if not lines[i+1]:
                print >> sys.stderr, "Error: Empty line:" + i + "\n"
                exit(1)
            if not lines[i+2].startswith("+"):
                print >> sys.stderr, "Error: Unusual line:" + i+2 + "\n", lines[i+2]
                exit(1)
            if not lines[i+3]:
                print >> sys.stderr, "Error: Empty line:" + 3 + "\n"
                exit(1)
            fastq.append([lines[i], lines[i+1], lines[i+2], lines[i+3]])
        if seed is not None:
            random.seed(seed)
        if seed:
            print seed
        downsample = random.sample(fastq, min(len(fastq),number))
        with gzip.open(outfile, "w") as f:
            for seq in downsample:
                for line in seq:
                    f.write(line+"\n")
    if os.path.isfile(outfile):
        return outfile


#### TOOLS #############################################################################################################

#### FASTQ DOWNSAMPLING ################################################################################################

def fastq_downsampling(args, indir):
    """
    Reduce the number of sequences in FASTQ file to n first sequences.
    """
    analysis_name = "FASTQ_downsampling"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        print "n:", args.fastq_downsample

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])

        logfile = "LOG"
        with open(logfile, "w") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        for infile in infiles:
            print "  {}".format(os.path.basename(infile))
            #print args.fastq_downsample
            #outfile = downsample_fastq(infile, os.path.basename(infile), number=args.fastq_downsample, rnd=True, seed=args.seed)     # true random downsampling
            outfile = downsample_fastq(infile, os.path.basename(infile), number=args.fastq_downsample)                               # downsampling from the head of the file
            if outfile:
                with open(logfile, "a+") as log:
                    log.write("{} downsampled to {} reads (from the beginning of the file)".format(os.path.basename(infile), args.fastq_downsample)+"\n")
            else:
                print "Error! Subsampling failed."
                exit(1)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### FASTQC ############################################################################################################

def run_fastqc(args, q, indir):
    """
    Run FastQC on FASTQ files.
    """
    analysis_name = "FastQC"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])

        logfile = "LOG"
        with open(logfile, "w") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        # FastQC
        for infile in infiles:

            jobs = ["{} --version".format(fastqc_path),
                    "{} --extract -t {} -o {} {}".format(fastqc_path, max(1,int(args.processors/args.threads)), os.getcwd(), infile)]

            q.put(Qjob(jobs, "LOG"))
            time.sleep(0.1)
        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### Trim Galore! ######################################################################################################

def run_trim_galore(args, q, indir):
    """
    Run Trim Galore! with user specified options.
    """
    analysis_name = "Trim Galore"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name.replace(" ", "_"))
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    #print "Outdir:", outdir
    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        logfile = "LOG"
        with open(logfile, "w") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        for infile in infiles:
            if args.paired:
                bname = re.sub("[1|2].fastq.gz$","",os.path.basename(infile[0]))
                jobs = ["{} --paired --fastqc {} {} {}".format(trim_galore_path, re.sub("[\"|\']$","",args.trim_galore_opts), infile[0], infile[1]),
                                "mv {}1_val_1.fq.gz {}1.fastq.gz".format(bname, bname),
                                "mv {}2_val_2.fq.gz {}2.fastq.gz".format(bname, bname)]
            else:
                bname = re.sub(".fastq.gz$","",os.path.basename(infile))
                jobs = ["{} --fastqc {} {}".format(trim_galore_path, args.trim_galore_opts, infile),
                        "mv {}_trimmed.fq.gz {}.fastq.gz".format(bname, bname)]

            q.put(Qjob(jobs, "LOG"))
            time.sleep(0.2)
        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### strand specificity and distance metrics ###########################################################################

def run_strand_specificity_and_distance_metrics(args, q, indir):
    """
    1) Downsampling to 1,000,000 reads
    2) Bowtie2 mapping
    5) infer_experiment (RSeQC) for strand specificity and inner_distance (RSeQC)
    6) Save a file with settings for TopHat2 (*.TopHat.txt): library-type, mate-inner-dist, mate-std-dev
    """
    analysis_name = "strand_specificity_and_distance_metrics"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        logfile = "LOG"
        with open(logfile, "a+") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        print "In:", os.path.abspath(indir)

        ## downsampling to 1000000 sequences
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")
        if args.paired:
            for pair in infiles:
                print " ".join(map(lambda x: os.path.basename(x),pair))
                for infile in pair:
                    bname = re.sub(".fastq.gz$","",os.path.basename(infile))
                    print "Downsample size:", downsample_size
                    print "Seed:", args.seed
                    downsample_fastq(infile, "{}.downsampled.fastq.gz".format(bname), number=downsample_size, rnd=True, seed=args.seed)
                print "FASTQ downsampled\n".format(downsample_size)
                with open(logfile, "a+") as log:
                    log.write(" ".join(map(lambda x: os.path.basename(x),pair))+"\n")
                    log.write("{} reads downsampled\n\n".format(downsample_size))
        else:
            for infile in infiles:
                print os.path.basename(infile)
                bname = re.sub(".fastq.gz$","",os.path.basename(infile))
                #downsample_fastq(infile, "{}.downsampled.fastq.gz".format(bname), number=downsample_size, rnd=True, seed=args.seed)
                downsample_fastq(infile, "{}.downsampled.fastq.gz".format(bname), number=downsample_size)
            print "FASTQ file(s) downsampled\n".format(downsample_size)
            with open(logfile, "a+") as log:
                #log.write("".join(map(lambda x: os.path.basename(x),infile))+"\n")
                log.write(infile+"\n")
                log.write("{} reads downsampled\n\n".format(downsample_size))

        ## Bowtie2 mapping
        infiles = check_for_paired_infiles(args, os.getcwd(), ".downsampled.fastq.gz")
        if args.paired:
            for pair in infiles:
                bname = re.sub("_R*[1|2].downsampled.fastq.gz$","",os.path.basename(pair[0]))
                cmd = "{} -x {} -1 {} -2 {} -p {} {} | {} view -Sb - | {} sort - {}.downsampled"\
                        .format(bowtie2_path, args.bowtie2_index, pair[0], pair[1], max(1,int(args.processors/args.threads)), args.bowtie_opts,
                        samtools_path,
                        samtools_path, bname)
                print cmd
                p = subprocess.check_call(cmd, shell=True)
                if p:
                    print "Error! Bowtie2 closed unexpectedly."
                    exit(1)
        else:
            for infile in infiles:
                bname = re.sub(".downsampled.fastq.gz$","",os.path.basename(infile))
                cmd = "{} -x {} -U {} -p {} {} | {} view -Sb - | {} sort - {}.downsampled"\
                        .format(bowtie2_path, args.bowtie2_index, infile, max(1,int(args.processors/args.threads)), args.bowtie_opts,
                        samtools_path,
                        samtools_path, bname)
                print cmd
                p = subprocess.check_call(cmd, shell=True)
                if p:
                    print "Error! Bowtie2 closed unexpectedly."
                    exit(1)

        infiles = sorted([f for f in os.listdir(os.getcwd()) if f.endswith(".downsampled.bam")])

        ## RSeQC infer_experiment (strand specificity)
        if not os.path.isdir("infer_experiment"):
            os.mkdir("infer_experiment")

            jobs = ["{} --version".format(os.path.join(rseqc_dir_path, "infer_experiment.py"))]
            q.put(Qjob(jobs, shell=False, logfile="LOG"))
            q.join()

            for infile in infiles:
                jobs = ["{} -i {} -r {} > infer_experiment/{}".format(os.path.join(rseqc_dir_path, "infer_experiment.py"), infile, args.bed, re.sub(".bam$","",os.path.basename(infile))+".infer_experiment.txt")]
                q.put(Qjob(jobs, shell=True, logfile="LOG"))
            time.sleep(0.1)
        q.join()

        ## Make settings file for TopHat
        for infile in infiles:
            bname = re.sub(".downsampled.bam$","",os.path.basename(infile))
            library_type = get_strand_from_rseqc("infer_experiment/{}.downsampled.infer_experiment.txt".format(bname), "tophat")
            with open("{}.TopHat2.txt".format(bname), "w") as f:
                f.write("library-type\t{}\n".format(library_type))

        if args.paired:
            ## Picard
            if args.insert_metrics == "Picard":

                ## Picard: CollectInsertSizeMetrics (mean, sd)
                if not os.path.isdir("InsertSizeMetrics"):
                    os.mkdir("InsertSizeMetrics")
                    for infile in infiles:
                        bname = re.sub(".downsampled.bam$","",os.path.basename(infile))

                        cmd = "java -jar -Xmx4g {}/CollectInsertSizeMetrics.jar INPUT={} OUTPUT=InsertSizeMetrics/{} HISTOGRAM_FILE=InsertSizeMetrics/{}".format(picardtools_dir_path, infile, bname+".InsertSizeMetrics.txt", bname+".histogram.pdf" )
                        print cmd
                        p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stderrdata = p.stderr.read()
                        print stderrdata
                        with open(logfile, "a+") as log:
                            log.write(cmd+"\n\n")
                            log.write(stderrdata+"\n")
                        p.wait()
                ## CollectAlignmentSummaryMetrics (mean read length)
                if not os.path.isdir("AlignmentSummaryMetrics"):
                    os.mkdir("AlignmentSummaryMetrics")
                    for infile in infiles:
                        bname = re.sub(".downsampled.bam$","",os.path.basename(infile))

                        cmd = "java -jar -Xmx4g {}/CollectAlignmentSummaryMetrics.jar INPUT={} OUTPUT=AlignmentSummaryMetrics/{}".format(picardtools_dir_path, infile, bname+".AlignmentSummaryMetrics.txt" )
                        print cmd
                        p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stderrdata = p.stderr.read()
                        print stderrdata
                        with open(logfile, "a+") as log:
                            log.write(cmd+"\n\n")
                            log.write(stderrdata+"\n")
                        p.wait()

                for infile in infiles:
                    bname = re.sub(".downsampled.bam$","",os.path.basename(infile))

                    with open("InsertSizeMetrics/{}.InsertSizeMetrics.txt".format(bname), "r") as f:
                        metrics = f.readlines()[7].split()
                        mean_insert_size = float(metrics[4])
                        standard_deviation = float(metrics[5])
                    with open("AlignmentSummaryMetrics/{}.AlignmentSummaryMetrics.txt".format(bname), "r") as f:
                        mean_read_length = float(f.readlines()[9].split()[15])
                    mate_inner_dist = "{:.0f}".format(mean_insert_size - mean_read_length*2)
                    mate_std_dev = "{:.0f}".format(standard_deviation)
                    with open("{}.TopHat2.txt".format(bname), "a") as f:
                        f.write("mate-inner-dist\t{}\n".format(mate_inner_dist))
                        f.write("mate-std-dev\t{}\n".format(mate_std_dev))

            ## RSeQC (default)
            else:
                ## RSeQC: inner_distance (mean, sd)
                if not os.path.isdir("inner_distance"):
                    os.mkdir("inner_distance")
                    for infile in infiles:
                        bname = re.sub(".bam$","",os.path.basename(infile))
                        jobs = ["{} -i {} -o inner_distance/{} -r {} > inner_distance/{}.inner_distance.summary.txt".format(os.path.join(rseqc_dir_path, "inner_distance.py"), infile, re.sub(".bam$","",os.path.basename(infile)), args.bed, bname) ]
                        q.put(Qjob(jobs, logfile="LOG", shell=True))
                q.join()

                for infile in infiles:
                    bname = re.sub(".downsampled.bam$","",os.path.basename(infile))
                    with open("inner_distance/{}.downsampled.inner_distance.summary.txt".format(bname), "r") as f:
                        lines = f.readlines()
                        mate_inner_dist = "{:.0f}".format(float(lines[1]))
                        mate_std_dev = "{:.0f}".format(float(lines[5]))
                    with open("{}.TopHat2.txt".format(bname), "a") as f:
                        f.write("mate-inner-dist\t{}\n".format(mate_inner_dist))
                        f.write("mate-std-dev\t{}\n".format(mate_std_dev))

        if is_error:
            exit(is_error)
        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### TOPHAT ############################################################################################################

def run_tophat(args, q, indir):
    """
    Run TopHat mapping.
    """
    analysis_name = "TopHat"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        logfile = "LOG"
        with open(logfile, "w") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        if args.paired:
            for pair in infiles:
                bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))

                #read metrics
                try:
                    tophat_settings_folder = [ d for d in os.listdir(args.outdir) if d.endswith("strand_specificity_and_distance_metrics") ]
                    if len([ d for d in os.listdir(args.outdir) if d.endswith("strand_specificity_and_distance_metrics") ]) != 1:
                        print "Error! Multiple folders that could contain TopHat settings files: {}".format(" ".join(tophat_settings_folder) )
                        exit(1)
                    tophat_settings_file = os.path.join(args.outdir, "{}/{}.TopHat2.txt".format(tophat_settings_folder[0], bname))
                    with open("{}".format(tophat_settings_file), "r") as f:
                        library_type = f.readline().split()[1]
                        mate_inner_dist = f.readline().split()[1]
                        mate_std_dev = f.readline().split()[1]
                except:
                    print "Error! Unable to read metrics from file: {}\n\n".format(tophat_settings_file)
                    exit(1)

                jobs = ["{} {} --num-threads {} --library-type {} --mate-inner-dist {} --mate-std-dev {} --output-dir {} --transcriptome-index {} {} {} {}".format(tophat2_path,
                                args.tophat_opts, max(1,int(args.processors/args.threads)), library_type, mate_inner_dist, mate_std_dev,
                                re.sub("_R*[1|2].fastq.gz","",os.path.basename(pair[0])), args.tophat2_index, args.bowtie2_index, pair[0], pair[1])]

                q.put(Qjob(jobs, "LOG"))
                time.sleep(0.1)
        else:
            for infile in infiles:
                bname = re.sub(".fastq.gz$", "", os.path.basename(infile))

                #read metrics
                try:
                    tophat_settings_folder = [ d for d in os.listdir(args.outdir) if d.endswith("strand_specificity_and_distance_metrics") ]
                    if len([ d for d in os.listdir(args.outdir) if d.endswith("strand_specificity_and_distance_metrics") ]) != 1:
                        print "Error! Multiple folders that could contain TopHat settings files: {}".format(" ".join(tophat_settings_folder) )
                        exit(1)
                    tophat_settings_file = os.path.join(args.outdir, "{}/{}.TopHat2.txt".format(tophat_settings_folder[0], bname))
                    with open("{}".format(tophat_settings_file), "r") as f:
                        library_type = f.readline().split()[1]
                except:
                    print "Error! Unable to read metrics from file: {}\n\n".format(tophat_settings_file)
                    exit(1)

                jobs = ["{} {} --num-threads {} --library-type {} --output-dir {} --transcriptome-index {} {} {}".format(tophat2_path,
                                args.tophat_opts,max(1,int(args.processors/args.threads)), library_type, bname, args.tophat2_index, args.bowtie2_index, infile)]

                q.put(Qjob(jobs, "LOG"))
                time.sleep(0.1)
        print
        q.join()
        if is_error:
            exit(is_error)

        ## Generate links in main TopHat output folder and index files
        for infile in infiles:
            if args.paired:
                bname = re.sub("_R*[1|2].fastq.gz$", "", os.path.basename(infile[0]))
            else:
                bname = re.sub(".fastq.gz$", "", os.path.basename(infile))
            tophat_file = "{}/accepted_hits.bam".format(bname)
            os.symlink(tophat_file, bname+".bam")
            subprocess.call("{} index {}.bam".format(samtools_path, bname), shell=True)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### RSeqQC ############################################################################################################

def run_rseqc(args, q, indir):
    """
    Run RSeQC with user specified options.
    """
    analysis_name = "RSeQC"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        print "In:", os.path.abspath(indir)

        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".bam")])

        logfile = "LOG"
        with open(logfile, "a+") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        if not os.path.isdir("infer_experiment"):
            os.mkdir("infer_experiment")
            for infile in infiles:
                jobs = ["{} -i {} -r {} > infer_experiment/{}".format(os.path.join(rseqc_dir_path, "infer_experiment.py"), infile, args.bed, re.sub(".bam$","",os.path.basename(infile))+".infer_experiment.txt")]
                q.put(Qjob(jobs, shell=True, logfile="LOG"))
            time.sleep(0.1)
            q.join()

        if not os.path.isdir("bam2wig"):
            os.mkdir("bam2wig")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))

                strand_rule = get_strand_from_rseqc("infer_experiment/{}.infer_experiment.txt".format(bname),"rseqc")
                if strand_rule is not None:
                    jobs = ["{} --strand='{}' -t 1000000000 --skip-multi-hits -i {} -o bam2wig/{} -s <(cut -f 1,2 {})".format(os.path.join(rseqc_dir_path, "bam2wig.py"), strand_rule, infile, bname, args.fasta_index) ]
                else:
                    jobs = ["{} -t 1000000000 --skip-multi-hits -i {} -o bam2wig/{} -s <(cut -f 1,2 {})".format(os.path.join(rseqc_dir_path, "bam2wig.py"), infile, bname, args.fasta_index) ]
                q.put(Qjob(jobs, logfile="LOG", shell=True))
            q.join()
            ## convert wig to bw
            for file in sorted([ os.path.join("bam2wig", f) for f in os.listdir("bam2wig") if f.endswith(".wig")]):
                bname = re.sub(".wig$","",os.path.basename(file))
                jobs = ["{} {} <( cut -f 1,2 {} ) bam2wig/{}.bw".format(os.path.join(ucsctools_dir_path, "wigToBigWig"), file, args.fasta_index, bname),
                        "rm {}".format(file)]
                q.put(Qjob(jobs, logfile="LOG", shell=True))
            q.join()

        if not os.path.isdir("bam_stat"):
            os.mkdir("bam_stat")
            for infile in infiles:
                jobs = ["{} -i {} 2> bam_stat/{}".format(os.path.join(rseqc_dir_path, "bam_stat.py"), infile, re.sub(".bam$","",os.path.basename(infile))+".bam_stat.txt")]
                q.put(Qjob(jobs, shell=True, logfile="LOG"))

        if not os.path.isdir("geneBody_coverage"):
            os.mkdir("geneBody_coverage")
            for infile in infiles:
                jobs = ["{} -i {} -o geneBody_coverage/{} -r {}".format(os.path.join(rseqc_dir_path, "geneBody_coverage.py"), infile, re.sub(".bam$","",os.path.basename(infile)), args.bed) ]
                q.put(Qjob(jobs, logfile="LOG"))

        if not os.path.isdir("junction_annotation"):
            os.mkdir("junction_annotation")
            for infile in infiles:
                jobs = ["{} -i {} -o junction_annotation/{} -r {}".format(os.path.join(rseqc_dir_path, "junction_annotation.py"), infile, re.sub(".bam$","",os.path.basename(infile)), args.bed) ]
                q.put(Qjob(jobs, logfile="LOG"))

        if not os.path.isdir("junction_saturation"):
            os.mkdir("junction_saturation")
            for infile in infiles:
                jobs = [ "{} -i {} -o junction_saturation/{} -r {}".format(os.path.join(rseqc_dir_path, "junction_saturation.py"), infile, re.sub(".bam$","",os.path.basename(infile)), args.bed) ]
                q.put(Qjob(jobs, logfile="LOG"))

        if not os.path.isdir("read_duplication"):
            os.mkdir("read_duplication")
            for infile in infiles:
                bname = re.sub(".bam","",os.path.basename(infile))
                jobs = [ "{} -i {} -o read_duplication/{}".format(os.path.join(rseqc_dir_path, "read_duplication.py"), infile, bname) ]
                q.put(Qjob(jobs, logfile="LOG", shell=True))

        ## WARNING this function consumes more than 8 GB RAM!!!
        if socket.gethostname() != "pc196":
            if not os.path.isdir("read_distribution"):
                os.mkdir("read_distribution")
                for infile in infiles:
                    jobs = [ "{} -i {} -r {} > read_distribution/{}.txt".format(os.path.join(rseqc_dir_path, "read_distribution.py"), infile, args.bed, re.sub(".bam$", "", os.path.basename(infile))) ]
                    q.put(Qjob(jobs, logfile="LOG", shell=True))

        if not os.path.isdir("read_GC"):
            os.mkdir("read_GC")
            for infile in infiles:
                jobs = [ "{} -i {} -o read_GC/{}".format(os.path.join(rseqc_dir_path, "read_GC.py"), infile, re.sub(".bam$", "", os.path.basename(infile))) ]
                q.put(Qjob(jobs, logfile="LOG"))

        if not os.path.isdir("read_NVC"):
            os.mkdir("read_NVC")
            for infile in infiles:
                jobs = [ "{} -i {} -o read_NVC/{}".format(os.path.join(rseqc_dir_path, "read_NVC.py"), infile, re.sub(".bam$", "", os.path.basename(infile))) ]
                q.put(Qjob(jobs, logfile="LOG"))

        ## DOES NOT WORK (for unknown reason; maybe because there are no quality values from TopHat? Same issue like in clipping_profile.py??)
        ##if not os.path.isdir("read_quality"):
        ##    os.mkdir("read_quality")
        ##q.put( "{} -i {} -o read_quality/{}".format(os.path.join(rseqc_dir_path, "read_quality.py"), infile, re.sub(".bam$", "", os.path.basename(infile))) ) #read_quality
        ## read_hexamer.py works only on fasta therefore skipped here; calculates hexamer frequencies

        if not os.path.isdir("spilt_bam"):
            os.mkdir("spilt_bam")
            for infile in infiles:
                jobs = [ "{} -i {} -r {} -o spilt_bam/{}".format(os.path.join(rseqc_dir_path, "split_bam.py"), infile, args.bed, re.sub(".bam$", "", os.path.basename(infile))) ]
                q.put(Qjob(jobs, logfile="LOG"))

        if args.paired:
            if not os.path.isdir("inner_distance"):
                os.mkdir("inner_distance")
                for infile in infiles:
                    bname = re.sub(".bam$","",os.path.basename(infile))
                    jobs = ["{} -i {} -o inner_distance/{} -r {} > inner_distance/{}.inner_distance.summary.txt".format(os.path.join(rseqc_dir_path, "inner_distance.py"), infile, re.sub(".bam$","",os.path.basename(infile)), args.bed, bname) ]
                    q.put(Qjob(jobs, logfile="LOG", shell=True))

        if not os.path.isdir("RPKM_count"):
            os.mkdir("RPKM_count")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))
                strand_rule = get_strand_from_rseqc("infer_experiment/{}.infer_experiment.txt".format(bname),"rseqc")
                if strand_rule is not None:
                    jobs = [ "{} --strand='{}' --skip-multi-hits -i {} -r {} -o RPKM_count/{}.RPKM".format(os.path.join(rseqc_dir_path, "RPKM_count.py"), strand_rule, infile, args.bed, re.sub(".bam$","",os.path.basename(infile))) ]
                else:
                    jobs = [ "{} --skip-multi-hits -i {} -r {} -o RPKM_count/{}.RPKM".format(os.path.join(rseqc_dir_path, "RPKM_count.py"), infile, args.bed, re.sub(".bam$","",os.path.basename(infile))) ]
                q.put(Qjob(jobs, logfile="LOG", shell=True))

        if not os.path.isdir("RPKM_saturation"):
            os.mkdir("RPKM_saturation")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))
                strand_rule = get_strand_from_rseqc("infer_experiment/{}.infer_experiment.txt".format(bname),"rseqc")
                if strand_rule is not None:
                    jobs =  [ "{} --strand='{}' -i {} -r {} -o RPKM_saturation/{}".format(os.path.join(rseqc_dir_path, "RPKM_saturation.py"), strand_rule, infile, args.bed, re.sub(".bam$","",os.path.basename(infile))) ]
                else:
                    jobs =  [ "{} -i {} -r {} -o RPKM_saturation/{}".format(os.path.join(rseqc_dir_path, "RPKM_saturation.py"), infile, args.bed, re.sub(".bam$","",os.path.basename(infile))) ]
                q.put(Qjob(jobs, logfile="LOG", shell=True))
        q.join()

        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### HTSEQ-COUNT #######################################################################################################

def run_htseq_count(args, q, indir):
    """
    Run htseq-count with user specified options.
    """
    analysis_name = "htseq-count"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".bam")])

        logfile = "LOG"
        with open(logfile, "w") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        for infile in infiles:
            bname = re.sub(".bam$","",os.path.basename(infile))
            strand_rule_folder = [ os.path.join(args.outdir,d,"infer_experiment") for d in os.listdir(args.outdir) if d.endswith("RSeQC") ][0]
            strand_rule = get_strand_from_rseqc(os.path.join(strand_rule_folder, bname+".infer_experiment.txt"), "htseq")
            jobs = [ "{samtools_path} sort -n {infile} -o {bam_outfile} -m {mem}G -@ {processors} | {samtools_path} view -h - | {htseq_count_path} {htseq_count_opts} --stranded={strand} - {gtf} > {outfile}" \
                      .format(samtools_path=samtools_path,
                              infile=infile,
                              bam_outfile=os.path.basename(os.path.dirname(infile))+".bam",
                              mem=max(10, int(samtools_max_mem/2)),
                              processors=2,
                              htseq_count_path=htseq_count_path,
                              strand=strand_rule,
                              htseq_count_opts=args.htseq_count_opts,
                              gtf=args.gtf,
                              outfile="{}.counts.txt".format(bname)) ]
            q.put(Qjob(jobs, logfile="LOG", shell=True))
            time.sleep(0.1)
        q.join()
        if is_error:
            exit(is_error)

        ## Combine count files into one file
        colnames = []
        count_dic = {}
        for infile in sorted([f for f in os.listdir(os.getcwd()) if f.endswith(".counts.txt")]):
            bname = re.sub(".counts.txt$","",os.path.basename(infile))
            colnames.append(bname)
            with open(infile, "r") as f:
                for line in f:
                    if not line.startswith("_"):
                        name, value = line.split()
                        if name not in count_dic:
                            count_dic[name] = [value]
                        else:
                            count_dic[name].append(value)
            with open("counts.txt", "w") as f:
                f.write("\t{}\n".format("\t".join(colnames)))
                for k,v in sorted(count_dic.iteritems()):
                    f.write("{}\t{}\n".format(k, "\t".join(v)))

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### featureCounts (Subread package) ###################################################################################

def run_featureCounts(args, q, indir):
    """
    Run featureCounts with user specified options.
    """
    analysis_name = "featureCounts"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".bam")])

        logfile = "LOG"
        with open(logfile, "w") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        jobs = ["{} -v".format(feature_counts_path)]
        q.put(Qjob(jobs, logfile="LOG", shell=False))
        q.join()

        for infile in infiles:
            bname = re.sub(".bam$","",os.path.basename(infile))
            print bname

            ## strand-specific counting
            strand_rule_folder = [ os.path.join(args.outdir,d,"infer_experiment") for d in os.listdir(args.outdir) if d.endswith("RSeQC") ][0]
            strand_rule = get_strand_from_rseqc(os.path.join(strand_rule_folder, bname+".infer_experiment.txt"), "htseq")
            if strand_rule == "reverse":
                strand_rule = 2
            elif strand_rule == "yes":
                strand_rule = 1
            else:
                strand_rule = 0

            #featureCounts -T 5 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.sam
            if args.paired:
                jobs = ["{} {} -p -B -C -T {} -s {} -a {} -o {} {}".format(feature_counts_path, args.featureCounts_opts, max(1,int(args.processors/args.threads)), strand_rule, args.gtf, "{}.counts.txt".format(bname), infile)]
                # -p : isPairedEnd
                # -B : requireBothEndsMap
                # -C : NOT countChimericFragments
                # -t : feature type; default: -t exon
            else:
                jobs = ["{} {} -C -T {} -a {} -o {} {}".format(feature_counts_path, args.featureCounts_opts, max(1,int(args.processors/args.threads)), args.gtf, "{}.counts.txt".format(bname), infile)]

            q.put(Qjob(jobs, logfile="LOG", shell=False))
            time.sleep(0.1)
        q.join()
        if is_error:
            exit(is_error)

        ## Combine count files into one file
        colnames = []
        count_dic = {}
        for infile in sorted([f for f in os.listdir(os.getcwd()) if f.endswith(".counts.txt")]):
            bname = re.sub(".counts.txt$","",os.path.basename(infile))
            colnames.append(bname)
            with open(infile, "r") as f:
                for line in f.readlines()[2:]:
                    elements = line.split()
                    name, value = elements[0], elements[6]
                    #print name, value
                    if name not in count_dic:
                        count_dic[name] = [value]
                    else:
                        count_dic[name].append(value)
            with open("counts.txt", "w") as f:
                f.write("\t{}\n".format("\t".join(colnames)))
                for k,v in sorted(count_dic.iteritems()):
                    f.write("{}\t{}\n".format(k, "\t".join(v)))

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### DESeq2 ############################################################################################################

def run_deseq2(args, q, indir):
    """
    Run DE analysis with DESeq2.
    """
    analysis_name = "DESeq2"
    args.analysis_counter += 1
    outdir = "{}_{}".format(add_leading_zeros(args.analysis_counter), analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)

        #print "In:", os.path.abspath(indir)
        #infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])
        counts_file = os.path.join(indir, "counts.txt")
        print "In:", counts_file

        logfile = "LOG"
        with open(logfile, "w") as log:
            log.write("Using {} threads for parallel file processing\n\n".format(args.threads))

        jobs = ["cat {} | /package/R-3.1.0/bin/R --vanilla --quiet --args {} {} {} {}".format(os.path.join(args.script,deseq2_path), args.sample_info, counts_file, 0.05, args.biomart)]

        q.put(Qjob(jobs, logfile="LOG", shell=True))

        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)



#### MAIN PROGRAM ######################################################################################################

def main():
    """
    Main program.
    """
    args = parse_args()
    ##print "Args:", args

    print "\n{} rna-seq-qc start".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'))
    if args.verbose:
        print "Genome:", args.genome
        print "Input dir:", args.indir
        print "Host:", socket.gethostname()
        print "Processors (max):", args.processors
        print "Threads (files in parallel):", args.threads
        print "Trim Galore options:", args.trim_galore_opts
        print "TopHat options:", args.tophat_opts
        print "htseq-count options:", args.htseq_count_opts
        print "Seed (random):", args.seed
        print "FASTA index:", args.fasta_index
        print "Bowtie2 index:", args.bowtie2_index
        print "Tophat2 index:", args.tophat2_index
        print "GTF:", args.gtf
        print "BED:", args.bed
        print "BioMart:", args.biomart

    ### Output dir
    ##print "Outdir:", args.outdir
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    os.chdir(args.outdir)

    args.analysis_counter = 0   # Counter for analyzes conducted
    indir = args.indir

    check_for_paired_infiles(args, indir, ".fastq.gz")      # sets the args.paired to True if sequences are paired end

    q = Queue()
    for i in range(args.threads):
        worker = Thread(target=queue_worker, args=(q, args.verbose))
        worker.setDaemon(True)
        worker.start()

    ## FASTQ downsampling
    if args.fastq_downsample:
        indir = fastq_downsampling(args, args.indir)
    #print "Indir:", indir

    ## Run FastQC
    run_fastqc(args, q, indir)

    # Run Trim Galore!
    if not args.no_trim:
        indir = run_trim_galore(args, q, indir)

    ## Stop here if requested by user (--no-bam)
    if args.no_bam:
        print "\nPipeline stopped because of user option '--no-bam'! rna-seq-qc finished (runtime: {})".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), datetime.datetime.now() - start)
        print "Output stored in: {}\n".format(args.outdir)
        exit(0)

    #print "Indir:", indir
    ## strand_specificity_and_distance_metrics
    run_strand_specificity_and_distance_metrics(args, q, indir)

    #print "Indir:", indir
    ## Run TopHat
    bam_dir = run_tophat(args, q, indir)

    ## Run RSeQC
    run_rseqc(args, q, bam_dir)

    ## Run htseq-count
    if args.count_prg == "featureCounts":
        count_dir = run_featureCounts(args, q, bam_dir)
    elif args.count_prg == "htseq-count":
        count_dir = run_htseq_count(args, q, bam_dir)
    else:
        print "Error! Unknown feature counting program:", args.count_prg
        exit(1)

    if args.sample_info:
        run_deseq2(args, q, count_dir)

    return args.outdir


if __name__ == "__main__":
    start = datetime.datetime.now()
    #print "os.environ['PATH']:    ", os.environ['PATH']
    #print "os.system('echo $PATH): ", os.system("echo $PATH")
    #print "Args:", sys.argv
    outdir = main()
    print "\n{} rna-seq-qc finished (runtime: {})".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), datetime.datetime.now() - start)
    print "Output stored in: {}\n".format(outdir)