#!/usr/bin/env python

__version__ = "rna-seq-qc v0.8.2"


__description__ = """
    {version}
    =================

    RNA-seq pipeline for processing RNA sequence data from high throughput
    sequencing.

    Fabian Kilpert - February 04, 2016
    email: kilpert@ie-freiburg.mpg.de

    This software is distributed WITHOUT ANY WARRANTY!
    --------------------------------------------------

    Following steps are executed in succession: FASTQ subsampling (optional),
    quality check with FASTQC, trimming of reads with Trim Galore (optional),
    estimation of insert size and strand specificity with RSeQC, mapping with
    TopHat, extensive quality check with RSeQC, counting of features with
    featureCounts (default) or htseq-count, differential expression analysis
    with DESeq2 (optional).

    The pipeline requires gzipped FASTQ files (.fastq.gz) for input, which are
    loaded from an input directory (-i INDIR). Read files belonging together
    require the exact same base name but ending either in "_R1" or "_R2" right
    in front of the .fastq.gz extension (e.g. reads_R1.fastq.gz,
    reads_R2.fastq.gz). In addition, a specific genome version argument must be
    provided (e.g. -g mm10) to define the reference data used for annotation.
    This loads a number of indexes for mapping programs (Bowtie2, TopHat2,
    HISAT2 (new!), etc.) from the corresponding configuration file of the
    rna-seq-qc sub-folder (e.g. rna-seq-qc/mm10.cfg). Additional genomes for
    selection can be provided as cfg-files by the user. The pipeline works for
    single end and paired end sequences alike.

    The DE analysis is only executed if a valid setup table is provided (e.g.
    --DE setup_table.tsv), which defines the relationships of the samples.

    More information on the pipeline can be found on the wiki page:
    http://epicenter/wiki/index.php/RNA-seq_pipeline_(rna-seq-pc.py)

    Example:
        python rna-seq-qc.py -i /path/to/fastq_dir -o /path/to/ouput_dir -g mm10 -v --DE sampleInfo.tsv
    """.format(version=__version__)

import argparse
from collections import OrderedDict
import datetime
# import gzip
import os
import os.path
from Queue import Queue
# import random
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import textwrap
from threading import Thread
import time


#### PATHS ####################################################################
temp_dir = tempfile.gettempdir()
script_path = os.path.dirname(os.path.realpath(__file__)) + "/"
python_path=sys.executable


#### Development settins ######################################################

if socket.gethostname() == "pc305.immunbio.mpg.de":

    ## Reading start options (if available)
    start_options_file = os.path.join(script_path, "rna-seq-qc", "start_options.txt")
    if os.path.isfile(start_options_file):
        start_options = []
        with open(start_options_file, "r") as f:
            for line in f.readlines():
                line = line.strip()
                if not line.startswith("#"):
                    if line:
                        start_options.extend(line.split())
        if start_options:
            print "Using start options:", " ".join(start_options)
            sys.argv = [sys.argv[0]] + start_options


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



## Path defaults to all needed scripts and programs + versions ################
fastqc_path = "/package/FastQC-0.11.3/"; fastqc_ver = "FastQC-0.11.3"
trim_galore_path = "/package/trim_galore_v0.4.0/"; trim_galore_ver = "TrimGalore-v0.4.0"
cutadapt_activate = "source /package/cutadapt-1.8.1/bin/activate &&"; cutadapt_ver = "Cutadapt-1.8.1"
rseqc_path = "/package/RSeQC-2.6.1/bin/"; rseqc_ver = "RSeQC-2.6.1"
rseqc_activate = "source /package/RSeQC-2.6.1/bin/activate &&"
bowtie2_path = "/package/bowtie2-2.2.3/"; bowtie2_ver = "Bowtie2-2.2.3"
bowtie2_export = "export PATH={}:$PATH &&".format(bowtie2_path)
picardtools_path = "/package/picard-tools-1.121/"; picardtools_ver = "Picard-tools-1.1.21"
tophat2_path = "/package/tophat-2.0.13.Linux_x86_64/"; tophat2_ver = "TopHat-2.0.13"
feature_counts_path = "/package/subread-1.5.0-p1/bin/"; feature_counts_ver = "featureCounts (subread-1.5.0-p1)"
htseq_count_path = "/package/HTSeq-0.6.1/bin/"; htseq_count_ver = "HTSeq-0.6.1"
R_path = "/package/R-3.2.0/bin/"; R_ver = "R-3.2.0"
samtools_path = "/package/samtools-1.2/"; samtools_ver = "Samtools-1.2"
samtools_export = "export PATH={}:$PATH &&".format(samtools_path)
ucsctools_dir_path = "/package/UCSCtools/"
#hisat_path = "/package/hisat-0.1.5-beta/bin/hisat"
#hisat_path = "/package/hisat-0.1.6-beta/bin/"; hisat_ver = "HISAT-0.1.6-beta"
hisat_path = "/package/hisat2-2.0.0-beta/"; hisat_ver = "HISAT2-2.0.0-beta"
R_libraries_export = "export R_LIBS_USER=/data/manke/repository/scripts/R/rna-seq-qc_libraries/R/x86_64-redhat-linux-gnu-library/3.2 &&"
deseq2_ver = "DESeq2-1.8.1"
deeptools_path = "/package/deeptools-2.0.0/bin/"; deeptools_ver = "deepTools-2.0"


## Different configurations for other physical machines
if socket.gethostname() == "pc305.immunbio.mpg.de":
    fastqc_path = "/home/kilpert/Software/bin/"; fastqc_ver = "FastQC"
    trim_galore_path = "/home/kilpert/Software/trim_galore/trim_galore_v0.4.0/"; trim_galore_ver = "TrimGalore-v0.4.0"
    cutadapt_activate = ""; cutadapt_ver = "Cutadapt"
    rseqc_path = "/home/kilpert/Software/RSeQC/RSeQC-2.6.3/build/scripts-2.7/"; rseqc_ver = "RSeQC-2.6.3"
    rseqc_activate = ""
    bowtie2_path = ""; bowtie2_ver = "Bowtie2"
    bowtie2_export = ""
    picardtools_path = "/home/kilpert/Software/picard-tools/picard-tools-1.115/"; picardtools_ver = "Picard-tools-1.115"
    tophat2_path = ""; tophat2_ver = "TopHat-2"
    feature_counts_path = ""; feature_counts_ver = "featureCounts"
    htseq_count_path = ""; htseq_count_ver = "HTSeq"
    R_path = "/usr/bin/"; R_ver = "R-3.2.3"
    samtools_path = ""; samtools_ver = "Samtools-1.2"
    samtools_export = ""
    ucsctools_dir_path = ""
    #hisat_path = "/home/kilpert/Software/hisat/hisat-0.1.5-beta/"; hisat_ver = "HISAT"
    hisat_path = "/home/kilpert/Software/hisat/hisat2-2.0.1-beta/"; hisat_ver = "HISAT2-2.0.1-beta"
    R_libraries_export = "export R_LIBS_USER=/data/manke/repository/scripts/rna-seq-qc/R/x86_64-pc-linux-gnu-library/3.2 &&"
    deseq2_ver = "DESeq2-1.8.1"
    deeptools_path = "/home/kilpert/Software/deepTools_release-1.6/deepTools/bin/"; deeptools_ver = "deepTools-1.6"


#### DEFAULT VARIABLES #################################################################################################
is_error = False
default_threads = 3           # Threads per process
default_parallel = 1          # Parallel files
samtools_mem = 1
samtools_threads = 1

if socket.gethostname().startswith("deep"):
    temp_dir = "/data/extended"
    default_threads = 4
    default_parallel = 6
    samtools_mem = 2
    samtools_threads = 2


#### COMMAND LINE PARSING ##############################################################################################

def parse_args():
    """
    Parse arguments from the command line.
    """
    parser = argparse.ArgumentParser(
    prog='rna-seq-qc.py',
    formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(__description__))

    ### Optional arguments
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument("-i", "--input-dir", dest="main_indir", required=True, help="Input dir containing (FASTQ)")
    parser.add_argument("-o", "--output-dir", dest="main_outdir", required=True, help="Output directory")
    parser.add_argument("--overwrite", dest="overwrite", action = "store_true", default=False, help="Overwrite results in existing folders!")
    parser.add_argument("-p", "--parallel", dest="parallel", metavar="INT", help="Number of files in parallel processing (default: {})".format(default_parallel), type=int, default=default_parallel)
    parser.add_argument("-t", "--threads", dest="threads", metavar="INT", help="Maximum number of threads for a single process (default: {})".format(default_threads), type=int, default=default_threads)
    parser.add_argument("--seed", dest="seed", metavar="INT", help="Random sampling seed", type=int, default=None)
    parser.add_argument("--fastq-downsample", dest="fastq_downsample", metavar="INT", help="Subsample first n fastq sequences (for testing only!)", type=int, default=None)
    parser.add_argument("-g", "--genome", dest="genome", required=True, help="Reference genome build")
    parser.add_argument("--trim_galore", dest="trim_galore_opts", metavar="STR", help="Trim Galore! option string (default: '--stringency 2')", type=str, default="--stringency 2")
    parser.add_argument("--tophat_opts", dest="tophat_opts", metavar="STR", help="TopHat2 option string", type=str, default="")     #--library-type fr-firststrand
    parser.add_argument("--bowtie_opts", dest="bowtie_opts", metavar="STR", help="Bowtie2 option string (default: '--end-to-end --fast'); for parameter estimation step only (no direct effect on Tophat2 mapping!)", type=str, default="--end-to-end --fast")
    parser.add_argument("--featureCounts_opts", dest="featureCounts_opts", metavar="STR", help="featureCounts option string (default: '')", type=str, default="-Q 10")
    parser.add_argument("--htseq-count_opts", dest="htseq_count_opts", metavar="STR", help="HTSeq htseq-count option string", type=str, default="--mode union")
    parser.add_argument("--insert-metrics", dest="insert_metrics", metavar="STR", help="Calculate insert size metrics (mean, sd) using Picard (default) or RSeQC", type=str, default="Picard")
    parser.add_argument("--count-prg", dest="count_prg", metavar="STR", help="Program used for counting features: featureCounts (default) or htseq-count", type=str, default="featureCounts")
    parser.add_argument("--library-type", dest="library_type", metavar="STR", help="Library type following TopHat naming scheme, e.g. fr-firstrand (default: auto)", type=str, default="auto")
    parser.add_argument("--mapping-prg", dest="mapping_prg", metavar="STR", help="Program used for mapping: TopHat2 (default) or HISAT2", type=str, default="TopHat2")
    parser.add_argument("--hisat_opts", dest="hisat_opts", metavar="STR", help="HISAT2 option string (default: '')", type=str, default="")
    parser.add_argument("--rseqc-preselection", dest="rseqc_preselection", help="Preselection of RSeQC programs; 1 (default) for minimum selection or 2 for maximum output", type=int, default="1")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose output")
    parser.add_argument("--report-secondary-alignments", dest="report_secondary_alignments", action="store_true", default=False, help="Output secondary alignments in BAM file. Default is to keep PRIMARY alignments only!!!")
    parser.add_argument("--bw", dest="bw", action="store_true", default=False, help="Generate BW (bigwig) files")
    parser.add_argument("--rseqc", dest="rseqc", action="store_true", default=False, help="Run RSeQC")
    parser.add_argument("--no-bam", dest="no_bam", action="store_true", default=False, help="First steps only. No alignment. No BAM file.")
    # parser.add_argument("--no-trim", dest="no_trim", action="store_true", default=False, help="Do not trim FASTQ reads. Default: Use Trim Galore! with default settings.")
    parser.add_argument("--trim", dest="trim", action="store_true", default=False, help="Activate trimming of fastq reads (default: no trimming)")
    parser.add_argument("--DE", dest="sample_info", help="Information on samples (required for DE analysis); see rna-seq-qc/example.setup_table.tsv for example. IMPORTANT: The first entry defines which group of samples are control. By this, the order of comparison and likewise the sign of values can be changed!")

    ### Positional arguments (required!)
    args = parser.parse_args()

    ## Tools path for report
    args.script_path = script_path

    ### Add useful paths
    args.cwd = os.getcwd()
    args.now = datetime.datetime.now()

    args.indir = args.main_indir
    args.outdir = args.main_outdir

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
    # if args.seed is None:
    #     args.seed = random.randint(12345, 987654321)

    ## Get reference data paths from config file
    ref_cfg_file_path = os.path.join(script_path, "rna-seq-qc/{}.cfg".format(args.genome))
    if not os.path.isfile(ref_cfg_file_path):
        print "Error! Configuration file NOT found for {}: {}".format(args.genome, ref_cfg_file_path)
        exit(1)
    configs = parse_config(ref_cfg_file_path)
    try:
        args.fasta_index = configs["fasta_index"]
        args.genome_index = configs["genome_index"]
        args.transcriptome_index = configs["transcriptome_index"]
        args.hisat_index = configs["hisat_index"]
        args.gtf = configs["gtf"]
        args.bed = configs["bed"]
    except:
        print "Error! Unable to read paths from config file:", ref_cfg_file_path
        exit(1)

    ## optional variables from config file
    if args.tophat_opts:
        #args.tophat_opts = re.sub( r'''^['"]+|['"]+$''', '', args.tophat_opts )     # strip quotes
        args.tophat_opts = args.tophat_opts.strip('"')
        args.tophat_opts = args.tophat_opts.strip("'")
    if args.hisat_opts:
        args.hisat_opts = args.hisat_opts.strip('"')
        args.hisat_opts = args.hisat_opts.strip("'")

    try:
        args.tophat_opts = args.tophat_opts + " " + configs["tophat_opts"]
        args.hisat_opts = args.hisat_opts + " " + configs["hisat_opts"]
    except:
        pass

    ## gene names, e.g. BioMart
    args.gene_names = os.path.join(script_path, "rna-seq-qc/{}.gene_names".format(args.genome))
    if not os.path.isfile(args.gene_names):
        args.gene_names = ""
        
    return args


#### GENERIC FUNCTIONS #################################################################################################

class Qjob():
    def __init__(self, cmds=None, cwd=None, logfile=None, backcopy=True, shell=False, keep_temp=False):
        self.cmds = cmds
        self.cwd = cwd
        self.logfile = logfile
        self.backcopy = backcopy
        self.shell = shell
        self.keep_temp = keep_temp
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
        m = re.match("^(.+)_R*[1|2]$", fname)
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


def run_subprocess(cmd, cwd, td, shell=False, logfile=None, backcopy=True, verbose=False, keep_temp=False):
    """
    Run the subprocess.
    """
    try:
        if verbose:
            print "Temp dir:", td, "\n"
            print cmd, "\n"

        sys.stdout.flush() # force printing

        if shell == True:
            return subprocess.check_call("cd {} && ".format(td)+cmd, shell=True, executable='/bin/bash')

        else:
            if type(cmd) is str:
                cmd = cmd.split()

            # print "CMD:", cmd
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=td)
            try:
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
                # print "poll:", p.poll()
                if p.poll() == None:
                    print "Error! Trying to terminate PID {}".format(p.pid)
                    p.terminate()
                    time.sleep(20)
                    if p.poll() == None:
                        print "Error! Trying to kill PID {}!".format(p.pid)
                        p.kill()
                    exit(1)
            return p.wait()
    finally:
        if backcopy:
            for f in os.listdir(td):
                if verbose:
                    print "Backcopy:", os.path.join(td, f)
                if os.path.isfile(os.path.join(cwd,f)) or os.path.isdir(os.path.join(cwd,f)):     # delete existing files with the same name in cwd
                    os.remove(os.path.join(cwd,f))
                #shutil.rmtree(os.path.join(cwd,f), ignore_errors=True)
                shutil.move(os.path.join(td,f), cwd)


def queue_worker(q, verbose=False, rest=0.2):
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
        td = tempfile.mkdtemp(prefix="rna-seq-qc.", dir=temp_dir)
        try:
            for cmd in job.cmds:
                if job.logfile:
                    with open(job.logfile, "a+") as log:
                        log.write("{}\n\n".format(cmd))
                return_code = run_subprocess(cmd, job.cwd, td, shell=job.shell, logfile=logfile, backcopy=job.backcopy, verbose=verbose, keep_temp=job.keep_temp)
                if return_code:
                    is_error = True
        finally:
            if os.path.isdir(td):
                if not job.keep_temp:
                    shutil.rmtree(td)
        q.task_done()
        time.sleep(rest)


def convert_library_type(library_type, prog, paired):
    """
    Convert TopHat library type to other program naming convention.
    """
    new = None

    if prog == "TopHat2":
        new == library_type


    elif prog == "RSeQC":
        if library_type == 'fr-firststrand' and paired==True:
            new = '"1+-,1-+,2++,2--"'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = '"1++,1--,2+-,2-+"'
        elif library_type == 'fr-firststrand' and paired==False:
            new = '"+-,-+"'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = '"++,--"'
        elif library_type == 'fr-unstranded':
            new = ''
        else:
            new = ''


    elif prog == "HISAT2":
        if library_type == 'fr-firststrand' and paired==True:
            new = 'RF'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = 'FR'
        elif library_type == 'fr-firststrand' and paired==False:
            new = 'R'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = 'F'
        elif library_type == 'fr-unstranded':
            new = 'unstranded'


    elif prog == "htseq-count":
        if library_type == 'fr-firststrand' and paired==True:
            new = 'reverse'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = 'yes'
        elif library_type == 'fr-firststrand' and paired==False:
            new = 'reverse'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = 'yes'
        elif library_type == 'fr-unstranded':
            new = 'no'


    elif prog == "featureCounts":
        if library_type == 'fr-firststrand' and paired==True:
            new = '2'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = '1'
        elif library_type == 'fr-firststrand' and paired==False:
            new = '2'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = '1'
        elif library_type == 'fr-unstranded':
            new = '0'

    return new



def get_library_type_from_rseqc(infile):
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

        threshold = 0.6      #min quotient threshold

        k = strands.keys()
        v = strands.values()

        specificity = "fr-unstranded"
        if '++,--' in strands.keys() and '+-,-+' in strands.keys():
            if strands['++,--'] >= threshold and strands['+-,-+'] <= threshold:
                specificity = "fr-secondstrand"
            elif strands['++,--'] <= threshold and strands['+-,-+'] >= threshold:
                specificity = "fr-firststrand"
        if '1++,1--,2+-,2-+' in strands.keys() and '1+-,1-+,2++,2--' in strands.keys():
            if strands['1++,1--,2+-,2-+'] >= threshold and strands['1+-,1-+,2++,2--'] <= threshold:
                specificity = "fr-secondstrand"
            elif strands['1++,1--,2+-,2-+'] <= threshold and strands['1+-,1-+,2++,2--'] >= threshold:
                specificity = "fr-firststrand"

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



def get_my_vars(infile, *var_names):
    """
    Read variables from file, return values.
    """
    values = {}
    with open(infile, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line:
                c = line.split()
                if len(c) != 2:
                    continue
                else:
                    if c[0] not in values:
                        values[c[0]] = c[1]
                    else:
                        print "Error! Variables names exist multiple times:", c[0]
                        exit(1)
    values_list = []
    for name in var_names:
        if name in values.keys():
            values_list.append(values[name])
        else:
            #print "Error! Requested variable does NOT exist in file:", name
            #exit(1)
            values_list.append(None)
    if len(values_list) == 1:
        return values_list[0]
    else:
        return values_list


#### TOOLS #############################################################################################################

#### FASTQ DOWNSAMPLING ################################################################################################

def run_fastq_downsampling(args, q, indir, analysis_name="FASTQ_downsampling"):
    """
    Reduce the number of sequences in FASTQ file to n first sequences
    """
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name.replace(" ", "_"))
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    print "Outdir:", outdir
    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])

        logfile = os.path.join(cwd, "LOG")
        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        for infile in infiles:
            if not args.seed:
                ##jobs = ["{} {}rna-seq-qc/downsample_fastq.py -v --head -n {} {} {}".format(python_path, script_path, args.fastq_downsample, infile, os.path.join(cwd, os.path.basename(infile)) ),]

                ## or:

                ## just from the head of the file using shell commands
                jobs = ["zcat {} 2>/dev/null | head -n{} | gzip > {}".format(infile, 4*int(args.fastq_downsample), os.path.join(cwd, os.path.basename(infile)) ),]
                ## Note that there is a not misleading "gzip: stdout: Broken pipe" message. The output is fine though!!!
            else:
                ## print "Using seed:", args.seed
                jobs = ["{} {}rna-seq-qc/downsample_fastq.py -v -s {} -n {} {} {}".format(python_path, script_path, args.seed, args.fastq_downsample, infile, os.path.join(cwd, os.path.basename(infile)) ),]

            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True) )

        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### FASTQC ############################################################################################################

def run_fastqc(args, q, indir, analysis_name="FastQC"):
    """
    Run FastQC on FASTQ files.
    """
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])

        logfile = os.path.join(cwd, "LOG")
        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        jobs = ["{}fastqc --version".format(fastqc_path)]
        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True))
        q.join()

        # FastQC
        for infile in infiles:
            jobs = ["{}fastqc --extract -o {} {}".format(fastqc_path, cwd, infile)]
            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True))
        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### Trim Galore! ######################################################################################################

def run_trim_galore(args, q, indir, analysis_name="Trim Galore"):
    """
    Run Trim Galore! with user specified options.
    """
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name.replace(" ", "_"))
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    #print "Outdir:", outdir
    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        logfile = os.path.join(cwd, "LOG")
        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        for infile in infiles:
            if args.paired:
                bname = re.sub("[1|2].fastq.gz$","",os.path.basename(infile[0]))
                jobs = ["{} {}trim_galore --paired {} {} {}".format(cutadapt_activate, trim_galore_path, args.trim_galore_opts, infile[0], infile[1]),
                                 "mv {}1_val_1.fq.gz {}1.fastq.gz".format(os.path.join(cwd,bname), os.path.join(cwd,bname)),
                                 "mv {}2_val_2.fq.gz {}2.fastq.gz".format(os.path.join(cwd,bname), os.path.join(cwd,bname))]
            else:
                bname = re.sub(".fastq.gz$","",os.path.basename(infile))
                jobs = ["{} {}trim_galore {} {}".format(cutadapt_activate, trim_galore_path, args.trim_galore_opts, infile),
                        "mv {}_trimmed.fq.gz {}.fastq.gz".format(os.path.join(cwd,bname), os.path.join(cwd,bname))]

            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True))
        q.join()
        if is_error:
            exit(is_error)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### library_type ###########################################################################

def run_library_type(args, q, indir):
    """
    - Random downsampling to n=500,000 reads
    - Bowtie2 mapping to genome -> library_type (infer_experiment; RSeQC)
    - Save a file with settings for TopHat2 (*.TopHat.txt): library-type
    """

    n = 500000  # number of downsampling reads

    analysis_name = "library_type"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        with open(logfile, "w") as log:
                log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        #######################################################################
        ## User defined library type
        #######################################################################
        if args.library_type != 'auto':
            library_type = args.library_type
            print "User defined library type:", library_type
            with open(logfile, "a") as log:
                log.write("User defined library type: {}\n\n".format(library_type))

        #######################################################################
        ## Automatic detection of library type
        #######################################################################
        else:
            print "Autodetecting library type..."

            ###################################################################
            ## downsampling
            ###################################################################

            if not args.seed:
                if args.paired:
                    for pair in infiles:
                        # print "PE"
                        bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))
                        jobs = ["bash {}rna-seq-qc/downsample_reservoir/downsample_reservoir_pe.sh {} {} {} {} {}".format(script_path,
                                    n,
                                    pair[0],
                                    pair[1],
                                    os.path.join(cwd,bname+"_R1.fastq.gz"),
                                    os.path.join(cwd,bname+"_R2.fastq.gz")),]
                        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, keep_temp=False))
                else:
                    for infile in infiles:
                        # print "SE"
                        bname = re.sub(".fastq.gz$","",os.path.basename(infile))
                        jobs = ["bash {}rna-seq-qc/downsample_reservoir/downsample_reservoir_se.sh {} {} {}".format(script_path,
                                    n,
                                    infile,
                                    os.path.join(cwd,bname+".fastq.gz")),]
                        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, keep_temp=False))
            else:
                infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])
                print infiles

                for infile in infiles:
                    jobs = ["{} {}rna-seq-qc/downsample_fastq.py -v -n {} -s {} {} {}".format(python_path, script_path, n, args.seed, infile, os.path.basename(infile) ),]

                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, keep_temp=False))

            q.join()
            if is_error:
                exit(is_error)


            ######################################################################
            ## Bowtie2 mapping to genome (for esimating strand specificity only!!!)
            ######################################################################
            print "CWD:", os.getcwd()
            infiles = check_for_paired_infiles(args, os.getcwd(), ".fastq.gz")

            if args.paired:
                for pair in infiles:
                    bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))

                    jobs = ["{}bowtie2 -x {} -1 {} -2 {} --threads {} {} | {}samtools view -Sb - | {}samtools sort -@ {} -m {}G - {}.genome_mapped"\
                             .format(bowtie2_path, args.genome_index, pair[0], pair[1], args.threads, args.bowtie_opts,
                             samtools_path,
                             samtools_path, samtools_threads, samtools_mem, bname),
                            "rm {} {}".format(pair[0], pair[1]),]

                    #q.put(Qjob(jobs, shell=True, logfile="LOG"))
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, shell=True))
            else:
                for infile in infiles:
                    bname = re.sub(".fastq.gz$","",os.path.basename(infile))

                    jobs = ["{}bowtie2 -x {} -U {} -p {} {} | {}samtools view -Sb - | {}samtools sort -@ {} -m {}G - {}.genome_mapped"\
                             .format(bowtie2_path, args.genome_index, infile, args.threads, args.bowtie_opts,
                             samtools_path,
                             samtools_path, samtools_threads, samtools_mem, bname),
                            "rm {}".format(infile),
                            ]

                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, shell=True))
            q.join()
            if is_error:
                exit(is_error)


            ###################################################################
            # RSeQC infer_experiment
            ###################################################################
            if not os.path.isdir("infer_experiment"):
                os.mkdir("infer_experiment")

            jobs = ["{} {}infer_experiment.py --version".format(rseqc_activate, rseqc_path)]
            q.put(Qjob(jobs, shell=True, logfile="LOG"))
            q.join()

            infiles = sorted([f for f in os.listdir(os.path.join(args.outdir, outdir)) if f.endswith(".genome_mapped.bam")])
            for infile in infiles:
                print infile
                bname = re.sub(".genome_mapped.bam$","",os.path.basename(infile))

                jobs = ["{} {}infer_experiment.py -i {} -r {} > {}"\
                            .format(rseqc_activate, rseqc_path, os.path.join(cwd, infile), args.bed, os.path.join(cwd, "infer_experiment", bname+".infer_experiment.txt")),
                        "rm {}".format(os.path.join(args.outdir, outdir, infile)),
                        ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, shell=True, keep_temp=False))
            q.join()

            ###################################################################
            ## combining results
            ###################################################################

            library_type_dic = {}

            infiles = sorted([f for f in os.listdir(os.path.join(args.outdir, outdir, "infer_experiment")) if f.endswith(".infer_experiment.txt")])
            for infile in infiles:
                bname = re.sub(".infer_experiment.txt$","",os.path.basename(infile))

                library_type = get_library_type_from_rseqc("infer_experiment/{}.infer_experiment.txt".format(bname))

                if library_type not in library_type_dic:
                    library_type_dic[library_type] = [bname]
                else:
                    library_type_dic[library_type].append(bname)

            print "Using TopHat2 naming scheme"

            for k,v in library_type_dic.iteritems():
                print "{} ({}): {}".format(k, len(v), " ".join(sorted(v)))

            library_type = sorted(library_type_dic.items(), key=lambda x: len(x[1]), reverse=True)[0][0]
            print "Auto library type:", library_type
            with open(logfile, "a") as log:
                log.write("Auto library type: {}\n\n".format(library_type))

        ## write librariy type(s) to file
        with open(os.path.join(cwd, "library_type.txt"), "w") as f:
            f.write("TopHat2\t{}\n".format(library_type))
            f.write("RSeQC\t{}\n".format(convert_library_type( library_type, 'RSeQC', args.paired ) ))
            f.write("HISAT2\t{}\n".format(convert_library_type( library_type, 'HISAT2', args.paired ) ))
            f.write("htseq-count\t{}\n".format(convert_library_type( library_type, 'htseq-count', args.paired ) ))
            f.write("featureCounts\t{}\n".format(convert_library_type( library_type, 'featureCounts', args.paired ) ))

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)



def run_distance_metrics(args, q, indir):
    """
    - Random downsampling to 1,000,000 reads for TopHat2
    - Bowtie2 mapping to transcriptome -> inner_distance (RSeQC) or InsertSizeMetrics (Picard) + CollectAlignmentSummaryMetrics (Picard)
    - Save a setting file: mate-inner-dist, mate-std-dev
    """
    n = 1000000  # number of downsampling reads

    analysis_name = "distance_metrics"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        with open(logfile, "w") as log:
                log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        ###################################################################
        ## downsampling
        ###################################################################
        if not args.seed:
            if args.paired:
                for pair in infiles:
                    # print "PE"
                    bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))
                    jobs = ["bash {}rna-seq-qc/downsample_reservoir/downsample_reservoir_pe.sh {} {} {} {} {}".format(script_path,
                                n,
                                pair[0],
                                pair[1],
                                os.path.join(cwd,bname+"_R1.fastq.gz"),
                                os.path.join(cwd,bname+"_R2.fastq.gz")),]
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, keep_temp=False))
            else:
                for infile in infiles:
                    # print "SE"
                    bname = re.sub(".fastq.gz$","",os.path.basename(infile))
                    jobs = ["bash {}rna-seq-qc/downsample_reservoir/downsample_reservoir_se.sh {} {} {}".format(script_path,
                                n,
                                infile,
                                os.path.join(cwd,bname+".fastq.gz")),]
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, keep_temp=False))
        else:
            infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".fastq.gz")])
            print infiles

            for infile in infiles:
                jobs = ["{} {}rna-seq-qc/downsample_fastq.py -v -n {} -s {} {} {}".format(python_path, script_path, n, args.seed, infile, os.path.basename(infile) ),]

                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, keep_temp=False))

        q.join()
        if is_error:
            exit(is_error)

        ######################################################################
        # Bowtie2 mapping to transcriptome
        ######################################################################

        infiles = check_for_paired_infiles(args, os.getcwd(), ".fastq.gz")

        ## Bowtie2 mapping to transcriptome (for insert size estimation!)
        for pair in infiles:
            bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))

            jobs = ["{}bowtie2 -x {} -1 {} -2 {} --threads {} {} | {}samtools view -Sb - | {}samtools sort -@ {} -m {}G - {}.transcriptome_mapped"\
                     .format(bowtie2_path, args.transcriptome_index, pair[0], pair[1], args.threads, args.bowtie_opts,
                     samtools_path,
                     samtools_path, samtools_threads, samtools_mem, bname),
                    "rm {} {}".format(pair[0], pair[1]),
                    ]

            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, backcopy=True, shell=True, keep_temp=False))
        q.join()
        if is_error:
            exit(is_error)

        #######################################################################

        infiles = sorted([f for f in os.listdir(os.getcwd()) if f.endswith(".transcriptome_mapped.bam")])

        ######################################################################
        ## Picard InsertSizeMetrics
        #######################################################################

        if args.insert_metrics == "Picard":

            ## Picard: CollectInsertSizeMetrics (mean, sd)
            if not os.path.isdir("InsertSizeMetrics"):
                os.mkdir("InsertSizeMetrics")

            for infile in infiles:
                bname = re.sub(".transcriptome_mapped.bam$","",os.path.basename(infile))

                jobs = ["export PATH={}:$PATH && java -jar -Xmx4g {}CollectInsertSizeMetrics.jar INPUT={} OUTPUT={} HISTOGRAM_FILE={}"\
                        .format( R_path, picardtools_path, os.path.join(cwd, infile),
                                 os.path.join(cwd, "InsertSizeMetrics", bname+".InsertSizeMetrics.txt"),
                                 os.path.join(cwd, "InsertSizeMetrics", bname+".histogram.pdf")),
                        ]

                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            if is_error:
                exit(is_error)

            ## CollectAlignmentSummaryMetrics (mean read length)
            if not os.path.isdir("AlignmentSummaryMetrics"):
                os.mkdir("AlignmentSummaryMetrics")

            for infile in infiles:
                bname = re.sub(".transcriptome_mapped.bam$","",os.path.basename(infile))

                jobs = ["java -jar -Xmx4g {}CollectAlignmentSummaryMetrics.jar INPUT={} OUTPUT={}"\
                            .format( picardtools_path,
                                     os.path.join(cwd, infile),
                                     os.path.join(cwd, "AlignmentSummaryMetrics", bname+".AlignmentSummaryMetrics.txt")),
                        #"rm {}".format(os.path.join(cwd, infile)),
                        ]

                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            if is_error:
                exit(is_error)

            for infile in infiles:
                bname = re.sub(".transcriptome_mapped.bam$","",os.path.basename(infile))

                with open("InsertSizeMetrics/{}.InsertSizeMetrics.txt".format(bname), "r") as f:
                    metrics = f.readlines()[7].split()
                    mean_insert_size = float(metrics[4])
                    standard_deviation = float(metrics[5])
                with open("AlignmentSummaryMetrics/{}.AlignmentSummaryMetrics.txt".format(bname), "r") as f:
                    mean_read_length = float(f.readlines()[9].split()[15])
                mate_inner_dist = "{:.0f}".format(mean_insert_size - mean_read_length*2)
                mate_std_dev = "{:.0f}".format(standard_deviation)
                with open("{}.TopHat2.txt".format(bname), "w") as f:
                    f.write("mate-inner-dist\t{}\n".format(mate_inner_dist))
                    f.write("mate-std-dev\t{}\n".format(mate_std_dev))
        else:
            ###############################################################
            ## RSeQC inner_distance
            ###############################################################
            if not os.path.isdir("inner_distance"):
                os.mkdir("inner_distance")

            ## RSeQC (default)
            if args.insert_metrics == "RSeQC":
                for infile in infiles:
                    bname = re.sub(".inner_distance.summary.txt$","",os.path.basename(infile))

                    jobs = ["export PATH={}:$PATH && {} {}inner_distance.py -i {} -o {} -r {} > {}".format(R_path, rseqc_activate, rseqc_path,
                                        os.path.join(cwd, infile),
                                        os.path.join(cwd, "inner_distance" , bname),
                                        args.bed,
                                        #bname) ]
                                        os.path.join(cwd, "inner_distance", bname+".inner_distance.summary.txt")),
                            #"rm {}".format(os.path.join(cwd, infile)),
                            ]
                    print jobs
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
                q.join()
                if is_error:
                    exit(is_error)

            for infile in infiles:
                bname = re.sub(".inner_distance.summary.txt$","",os.path.basename(infile))
                with open("inner_distance/{}.inner_distance.summary.txt".format(bname), "r") as f:
                    lines = f.readlines()
                    e = lines[1].split()
                    mate_inner_dist = "{:.0f}".format(float(e[1]))
                    mate_std_dev = "{:.0f}".format(float(e[3]))
                with open("{}.TopHat2.txt".format(bname), "w") as f:
                    f.write("mate-inner-dist\t{}\n".format(mate_inner_dist))
                    f.write("mate-std-dev\t{}\n".format(mate_std_dev))

    print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


#### TOPHAT ############################################################################################################

def run_tophat(args, q, indir):
    """
    Run TopHat mapping.
    """
    analysis_name = "TopHat2"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        ## read library type from file
        library_type = get_my_vars( os.path.join(args.outdir, "library_type", "library_type.txt"), "TopHat2" )

        ## set variable to filter for multi mapping reads
        if args.report_secondary_alignments:
            report_secondary_alignments = "--report-secondary-alignments"   # secondary alignments
        else:
            report_secondary_alignments = ""   # no secondary alignments (default)

        if args.paired:
            for pair in infiles:
                bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))

                ## read metrics from file
                mate_inner_dist, mate_std_dev = get_my_vars( os.path.join(args.outdir, "distance_metrics", bname+".TopHat2.txt"), "mate-inner-dist", "mate-std-dev" )

                jobs = ["{} {} {}tophat2 {} --num-threads {} {} --library-type {} --mate-inner-dist {} --mate-std-dev {} --output-dir {} --transcriptome-index {} {} {} {}"\
                            .format(bowtie2_export, samtools_export, tophat2_path, args.tophat_opts, args.threads, report_secondary_alignments, library_type, mate_inner_dist, mate_std_dev,
                            os.path.join(cwd, bname), args.transcriptome_index, args.genome_index, pair[0], pair[1])]
                            #re.sub("_R*[1|2].fastq.gz","",os.path.basename(pair[0])), args.transcriptome_index, args.genome_index, pair[0], pair[1])]

                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
        else:
            for infile in infiles:
                bname = re.sub(".fastq.gz$", "", os.path.basename(infile))

                jobs = ["{} {} {}tophat2 {} --num-threads {} {} --library-type {} --output-dir {} --transcriptome-index {} {} {}"\
                            .format(bowtie2_export, samtools_export, tophat2_path, args.tophat_opts, args.threads, report_secondary_alignments, library_type, os.path.join(cwd, bname), args.transcriptome_index, args.genome_index, infile)]

                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
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
            subprocess.call("{}samtools index {}.bam".format(samtools_path, bname), shell=True)

        print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)


def run_hisat2(args, q, indir):
    """
    Run HISAT2 mapping.
    """
    analysis_name = "HISAT2"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        print "In:", os.path.abspath(indir)
        infiles = check_for_paired_infiles(args, indir, ".fastq.gz")

        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        ## read library type from file
        library_type = get_my_vars( os.path.join(args.outdir, "library_type", "library_type.txt"), "HISAT2" )
        if library_type == "unstranded":
            library_type = ""
        else:
            library_type = "--rna-strandness "+library_type

        ## set variable to filter for multi mapping reads
        if args.report_secondary_alignments:
            report_secondary_alignments = ""         # all alignments
        else:
            report_secondary_alignments = "-F256"    # primary alignments only (default)

        ## PE
        if args.paired:
            for pair in infiles:
                bname = re.sub("_R*[1|2].fastq.gz$","",os.path.basename(pair[0]))

                if not os.path.isdir( os.path.join(cwd, bname) ):
                    os.mkdir( os.path.join(cwd, bname) )

                cmdl = "{}hisat2 {} -p {} -x {} {} -1 {} -2 {} --novel-splicesite-outfile {} --un-conc-gz {} --al-conc-gz {} --met-file {} 2> {} | {}samtools view {} -Sb - | {}samtools sort -@ {} -m {}G - {}"\
                            .format(hisat_path, args.hisat_opts, args.threads, args.hisat_index, library_type, pair[0], pair[1],
                                    os.path.join(cwd, bname+"/"+"splice_sites.txt"),
                                    os.path.join(cwd, bname+"/"+"un-conc.fastq.gz"),        # --un-conc
                                    os.path.join(cwd, bname+"/"+"al-conc.fastq.gz"),        # --al-conc
                                    os.path.join(cwd, bname+"/"+"metrics.txt"),
                                    os.path.join(cwd, bname+"/"+"align_summary.txt"),
                                    samtools_path,
                                    report_secondary_alignments,
                                    samtools_path, samtools_threads, samtools_mem,
                                    os.path.join(cwd, bname+"/"+"accepted_hits"),
                                    )

                jobs = [cmdl,
                        ## add command line to bam header
                        "cat <({}samtools view -H {}) <(echo '@PG\tCL:{}') > {}"\
                            .format(samtools_path, os.path.join(cwd, bname+"/"+"accepted_hits.bam"), cmdl, os.path.join(cwd, bname+"/"+"header.sam")),
                        "{}samtools reheader {} {} | {}samtools view -bS - > {}"\
                            .format(samtools_path, os.path.join(cwd, bname+"/"+"header.sam"), os.path.join(cwd, bname+"/"+"accepted_hits.bam"), samtools_path, os.path.join(cwd, bname+"/"+"accepted_hits.reheader.bam") ),
                        "mv {} {}".format( os.path.join(cwd, bname+"/"+"accepted_hits.reheader.bam"), os.path.join(cwd, bname+"/"+"accepted_hits.bam") ),
                        "rm {}".format(os.path.join(cwd, bname+"/"+"header.sam"), ),
                        ]

                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
        ## SE
        else:
            for infile in infiles:
                bname = re.sub(".fastq.gz$", "", os.path.basename(infile))

                if not os.path.isdir( os.path.join(cwd, bname) ):
                    os.mkdir( os.path.join(cwd, bname) )

                cmdl = "{}hisat2 {} -p {} -x {} {} -U {} --novel-splicesite-outfile {} --un-gz {} --al-gz {} --met-file {} 2> {} | {}samtools view {} -Sb - | {}samtools sort -@ {} -m {}G - {}"\
                            .format(hisat_path, args.hisat_opts, args.threads, args.hisat_index, library_type, infile,
                                    os.path.join(cwd, bname+"/"+"splice_sites.txt"),
                                    os.path.join(cwd, bname+"/"+"un.fastq.gz"),         # --un
                                    os.path.join(cwd, bname+"/"+"al.fastq.gz"),         # --al
                                    os.path.join(cwd, bname+"/"+"metrics.txt"),
                                    os.path.join(cwd, bname+"/"+"align_summary.txt"),
                                    samtools_path,
                                    report_secondary_alignments,
                                    samtools_path, samtools_threads, samtools_mem,
                                    os.path.join(cwd, bname+"/"+"accepted_hits"),
                                    )

                jobs = [cmdl,
                        # add command line to bam header
                        "cat <({}samtools view -H {}) <(echo '@PG\tCL:{}') > {}"\
                            .format(samtools_path, os.path.join(cwd, bname+"/"+"accepted_hits.bam"), cmdl, os.path.join(cwd, bname+"/"+"header.sam")),
                        "{}samtools reheader {} {} | {}samtools view -bS - > {}"\
                            .format(samtools_path, os.path.join(cwd, bname+"/"+"header.sam"), os.path.join(cwd, bname+"/"+"accepted_hits.bam"), samtools_path, os.path.join(cwd, bname+"/"+"accepted_hits2.bam") ),
                        "mv {} {}".format( os.path.join(cwd, bname+"/"+"accepted_hits2.bam"), os.path.join(cwd, bname+"/"+"accepted_hits.bam") ),
                        "rm {}".format(os.path.join(cwd, bname+"/"+"header.sam")),
                        ]

                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
        print
        q.join()
        if is_error:
            exit(is_error)

        ## Generate links in main HISAT2 output folder and index files
        for infile in infiles:
            if args.paired:
                bname = re.sub("_R*[1|2].fastq.gz$", "", os.path.basename(infile[0]))
            else:
                bname = re.sub(".fastq.gz$", "", os.path.basename(infile))
            tophat_file = "{}/accepted_hits.bam".format(bname)
            os.symlink(tophat_file, bname+".bam")
            subprocess.call("{}samtools index {}.bam".format(samtools_path, bname), shell=True)

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
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".bam")])

        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        ## read library type
        library_type = get_my_vars( os.path.join(args.outdir, "library_type", "library_type.txt"), "htseq-count" )

        for infile in infiles:
            bname = re.sub(".bam$","",os.path.basename(infile))

            jobs = [ "{samtools_path}samtools sort -n {infile} -o {bam_outfile} -@ {threads} -m {mem}G | {samtools_path}samtools view -h - | {htseq_count_path}htseq-count {htseq_count_opts} --stranded={strand} - {gtf} > {outfile}" \
                      .format(samtools_path=samtools_path,
                              infile=infile,
                              bam_outfile=os.path.basename(os.path.dirname(infile))+".bam",
                              threads=samtools_threads,
                              mem=samtools_mem,
                              htseq_count_path=htseq_count_path,
                              strand=library_type,
                              htseq_count_opts=args.htseq_count_opts,
                              gtf=args.gtf,
                              outfile="{}.counts.txt".format(bname)) ]
            #q.put(Qjob(jobs, logfile="LOG", shell=True))
            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
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
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".bam")])

        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        jobs = ["{}featureCounts -v".format(feature_counts_path)]
        q.put(Qjob(jobs, logfile="LOG", shell=False))
        q.join()


        ## read library type
        library_type = get_my_vars( os.path.join(args.outdir, "library_type", "library_type.txt"), "featureCounts" )

        for infile in infiles:
            bname = re.sub(".bam$","",os.path.basename(infile))

            #featureCounts -T 5 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.sam
            if args.paired:
                jobs = ["{}featureCounts {} -p -B -C -T {} -s {} -a {} -o {} {}".format(feature_counts_path, args.featureCounts_opts, args.threads, library_type, args.gtf, "{}.counts.txt".format(bname), infile)]
                # -p : isPairedEnd
                # -B : requireBothEndsMap
                # -C : NOT countChimericFragments
                # -t : feature type; default: -t exon
            else:
                jobs = ["{}featureCounts {} -C -T {} -a {} -o {} {}".format(feature_counts_path, args.featureCounts_opts, args.threads, args.gtf, "{}.counts.txt".format(bname), infile)]

            #q.put(Qjob(jobs, logfile="LOG", shell=False))
            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=False, backcopy=True, keep_temp=False))
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


#### BW files ############################################################################################################
def run_bigWig(args, q, indir):
    """
    Run generating BW files
    """
    analysis_name = "bigWig"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        print "In:", os.path.abspath(indir)
        infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".bam")])

        with open(logfile, "w") as log:
            log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

        ## PE
        if (args.paired and args.library_type == 'fr-firststrand'):
            for infile in infiles:
                ##bname = re.sub(".bam$","",os.path.basename(infile))

                bamCoverage_path = os.path.join(deeptools_path, "bamCoverage")
                if not os.path.isfile(bamCoverage_path):
                    print "Error! bamCoverage NOT found:"
                    print bamCoverage_path
                    exit(2)
                else:
                    jobs = ["bash {}rna-seq-qc/bw/RNA_bw_PE_best.sh {} {} {} {} {}".format(script_path, infile, os.path.join(args.outdir,outdir,"best"), deeptools_path, samtools_path, args.threads),
                            "bash {}rna-seq-qc/bw/RNA_bw_PE_all.sh {} {} {} {} {}".format(script_path, infile, os.path.join(args.outdir,outdir,"all"), deeptools_path, samtools_path, args.threads),]
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            if is_error:
                exit(is_error)
        ## SE
        else:
            ##print "Only available for PE + fr-firstrand data!"
            for infile in infiles:
                ##bname = re.sub(".bam$","",os.path.basename(infile))

                bamCoverage_path = os.path.join(deeptools_path, "bamCoverage")
                if not os.path.isfile(bamCoverage_path):
                    print "Error! bamCoverage NOT found:"
                    print bamCoverage_path
                    exit(2)
                else:
                    jobs = ["bash {}rna-seq-qc/bw/RNA_bw_SE_best.sh {} {} {} {} {}".format(script_path, infile, os.path.join(args.outdir,outdir,"best"), deeptools_path, samtools_path, args.threads),
                            "bash {}rna-seq-qc/bw/RNA_bw_SE_all.sh {} {} {} {} {}".format(script_path, infile, os.path.join(args.outdir,outdir,"all"), deeptools_path, samtools_path, args.threads),]
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            if is_error:
                exit(is_error)

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
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        counts_file = os.path.join(indir, "counts.txt")
        print "In:", counts_file

        jobs = ["{} cat {}rna-seq-qc/DESeq2.R | {}R --vanilla --quiet --args {} {} {} {}".format(R_libraries_export, script_path, R_path, args.sample_info, counts_file, 0.05, args.gene_names),
                "[ -f {rplotspath} ] && rm {rplotspath}".format( rplotspath=os.path.join(cwd,"Rplots.pdf") ),
                ]
        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))

        q.join()
        if is_error:
            exit(is_error)

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
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    os.chdir(outdir)
    cwd = os.getcwd()
    logfile = os.path.join(cwd, "LOG")

    #print "In:", os.path.abspath(indir)

    strand_rule = convert_library_type( get_my_vars(os.path.join(args.outdir,"library_type","library_type.txt"),"RSeQC"), 'RSeQC', args.paired )

    infiles = sorted([os.path.join(indir, f) for f in os.listdir(os.path.abspath(indir)) if f.endswith(".bam")])

    with open(logfile, "a+") as log:
        log.write("Processing {} file(s) in parallel\n\n".format(args.parallel))

    #print "RSeQC preselection:", args.rseqc_preselection

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("infer_experiment"):
            t1 = datetime.datetime.now()
            os.mkdir("infer_experiment")
            for infile in infiles:
                bname = ".".join(os.path.basename(infile).split(".")[:-1])
                jobs = ["export PATH={}:$PATH && {} {}infer_experiment.py -i {} -r {} > {}".format(R_path, rseqc_activate, rseqc_path, os.path.join(cwd, infile), args.bed, os.path.join(cwd, "infer_experiment", bname)+".infer_experiment.txt")]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 2:
        if not os.path.isdir("bam2wig"):
            t1 = datetime.datetime.now()
            os.mkdir("bam2wig")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))

                if strand_rule is not None:
                    jobs = ["{} {}bam2wig.py --strand='{}' -t 1000000000 --skip-multi-hits -i {} -o {} -s {}".format(rseqc_activate, rseqc_path, strand_rule, infile, os.path.join(cwd, "bam2wig", bname), args.fasta_index) ]
                else:
                    jobs = ["{} {}bam2wig.py -t 1000000000 --skip-multi-hits -i {} -o {} -s {}".format(rseqc_activate, rseqc_path, infile, os.path.join(cwd, "bam2wig", bname), args.fasta_index) ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            ## convert wig to bw
            for file in sorted([ os.path.join("bam2wig", f) for f in os.listdir("bam2wig") if f.endswith(".wig")]):
                bname = re.sub(".wig$","",os.path.basename(file))
                jobs = ["{}wigToBigWig {} <( cut -f 1,2 {} ) {}".format(ucsctools_dir_path, os.path.join(cwd, file), args.fasta_index, os.path.join(cwd, "bam2wig", bname)+".bw"),
                        "rm {}".format(os.path.join(cwd, file))]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("bam_stat"):
            t1 = datetime.datetime.now()
            os.mkdir("bam_stat")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))
                jobs = ["export PATH={}:$PATH && {} {}bam_stat.py -i {} 2> {}".format(R_path, rseqc_activate, rseqc_path, infile, os.path.join(cwd, "bam_stat", bname)+".bam_stat.txt")]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 2:
        if not os.path.isdir("geneBody_coverage"):
            t1 = datetime.datetime.now()
            os.mkdir("geneBody_coverage")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))
                jobs = ["export PATH={}:$PATH && {} {}geneBody_coverage.py -i {} -o {} -r {}".format(R_path, rseqc_activate, rseqc_path, infile, os.path.join(cwd, "geneBody_coverage", bname)+".geneBody_coverage", args.bed) ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("junction_annotation"):
            t1 = datetime.datetime.now()
            os.mkdir("junction_annotation")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))
                jobs = ["export PATH={}:$PATH && {} {}junction_annotation.py -i {} -o {} -r {}".format(R_path, rseqc_activate, rseqc_path, infile, os.path.join(cwd, "junction_annotation", bname)+".junction_annotation", args.bed) ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("junction_saturation"):
            t1 = datetime.datetime.now()
            os.mkdir("junction_saturation")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))
                jobs = [ "export PATH={}:$PATH && {} {}junction_saturation.py -i {} -o {} -r {}".format(R_path, rseqc_activate, rseqc_path, infile, os.path.join(cwd, "junction_saturation", bname)+".junction_saturation", args.bed) ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("read_duplication"):
            t1 = datetime.datetime.now()
            os.mkdir("read_duplication")
            for infile in infiles:
                bname = re.sub(".bam","",os.path.basename(infile))
                jobs = [ "export PATH={}:$PATH && {} {}read_duplication.py -i {} -o {}".format(R_path, rseqc_activate, rseqc_path, infile, os.path.join(cwd, "read_duplication", bname)+".read_duplication") ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    ## WARNING this function consumes more than 8 GB RAM!!!
    if args.rseqc_preselection >= 2:
        if socket.gethostname().startswith("deep"):
            if not os.path.isdir("read_distribution"):
                t1 = datetime.datetime.now()
                os.mkdir("read_distribution")
                for infile in infiles:
                    bname = re.sub(".bam","",os.path.basename(infile))
                    jobs = [ "export PATH={}:$PATH && {} {}read_distribution.py -i {} -r {} > {}".format(R_path, rseqc_activate, rseqc_path, infile, args.bed, os.path.join(cwd, "read_distribution", bname)+".read_distribution.txt") ]
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
                q.join()
                t2 = datetime.datetime.now()
                print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("read_GC"):
            t1 = datetime.datetime.now()
            os.mkdir("read_GC")
            for infile in infiles:
                bname = re.sub(".bam","",os.path.basename(infile))
                jobs = [ "export PATH={}:$PATH && {} {}read_GC.py -i {} -o {}".format(R_path, rseqc_activate, rseqc_path, infile, os.path.join(cwd, "read_GC", bname)+".read_GC") ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("read_NVC"):
            t1 = datetime.datetime.now()
            os.mkdir("read_NVC")
            for infile in infiles:
                bname = re.sub(".bam","",os.path.basename(infile))
                jobs = [ "export PATH={}:$PATH && {} {}read_NVC.py -i {} -o {}".format(R_path, rseqc_activate, rseqc_path, infile, re.sub(".bam$", "", os.path.join(cwd, "read_NVC", bname)+".read_NVC")) ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    ## DOES NOT WORK (for unknown reason; maybe because there are no quality values from TopHat? Same issue like in clipping_profile.py??)
    ##if not os.path.isdir("read_quality"):
    ##    os.mkdir("read_quality")
    ##q.put( "{}read_quality.py -i {} -o read_quality/{}".format(rseqc_path, infile, re.sub(".bam$", "", os.path.basename(infile))) ) #read_quality
    ## read_hexamer.py works only on fasta therefore skipped here; calculates hexamer frequencies

    if args.rseqc_preselection >= 2:
        if not os.path.isdir("spilt_bam"):
            t1 = datetime.datetime.now()
            os.mkdir("spilt_bam")
            for infile in infiles:
                bname = re.sub(".bam","",os.path.basename(infile))
                jobs = [ "export PATH={}:$PATH && {} {}split_bam.py -i {} -r {} -o {}".format(R_path, rseqc_activate, rseqc_path, infile, args.bed, re.sub(".bam$", "", os.path.join(cwd, "spilt_bam", bname)+".spilt_bam")) ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if args.paired:
            if not os.path.isdir("inner_distance"):
                t1 = datetime.datetime.now()
                os.mkdir("inner_distance")
                for infile in infiles:
                    bname = re.sub(".bam$","",os.path.basename(infile))
                    jobs = ["export PATH={}:$PATH && {} {}inner_distance.py -i {} -o {} -r {} > {}".format(R_path, rseqc_activate, rseqc_path, os.path.join(cwd, infile), os.path.join(cwd, "inner_distance", bname), args.bed, os.path.join(cwd, "inner_distance", bname)+".inner_distance.summary.txt") ]
                    q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
                q.join()
                t2 = datetime.datetime.now()
                print "Duration:", t2-t1

    if args.rseqc_preselection >= 1:
        if not os.path.isdir("RPKM_count"):
            t1 = datetime.datetime.now()
            os.mkdir("RPKM_count")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))

                if strand_rule:
                    jobs = [ "export PATH={}:$PATH && {} {}RPKM_count.py --strand='{}' --skip-multi-hits -i {} -r {} -o {}".format(R_path, rseqc_activate, rseqc_path, strand_rule, os.path.join(cwd, infile), args.bed,  os.path.join(cwd, "RPKM_count", bname)+ ".RPKM") ]
                else:
                    jobs = [ "export PATH={}:$PATH && {} {}RPKM_count.py --skip-multi-hits -i {} -r {} -o {}".format(R_path, rseqc_activate, rseqc_path, os.path.join(cwd, infile), args.bed, os.path.join(cwd, "RPKM_count", bname)+".RPKM") ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if args.rseqc_preselection >= 2:
        if not os.path.isdir("RPKM_saturation"):
            t1 = datetime.datetime.now()
            os.mkdir("RPKM_saturation")
            for infile in infiles:
                bname = re.sub(".bam$","",os.path.basename(infile))

                if strand_rule:
                    jobs =  [ "export PATH={}:$PATH && {} {}RPKM_saturation.py --strand='{}' -i {} -r {} -o {}".format(R_path, rseqc_activate, rseqc_path, strand_rule, os.path.join(cwd, infile), args.bed, os.path.join(cwd, "RPKM_saturation", bname)) ]
                else:
                    jobs =  [ "export PATH={}:$PATH && {} {}RPKM_saturation.py -i {} -r {} -o {}".format(R_path, rseqc_activate, rseqc_path, os.path.join(cwd, infile), args.bed, os.path.join(cwd, "RPKM_saturation", bname)) ]
                q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            t2 = datetime.datetime.now()
            print "Duration:", t2-t1

    if is_error:
        exit(is_error)

    #print "Out:", os.path.join(args.outdir, outdir)
    os.chdir(args.outdir)
    return os.path.join(args.outdir, outdir)



def run_project_report(args, q):
    """
    Compile a report on the Project.
    """
    analysis_name = "project_report"
    args.analysis_counter += 1
    outdir = "{}".format(analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)

    if args.overwrite and os.path.isdir(outdir):
        shutil.rmtree(outdir)

    if os.path.isdir(outdir):
        print "Output folder already present: {}".format(outdir)
    else:
        os.mkdir(outdir)
        os.chdir(outdir)
        cwd = os.getcwd()
        logfile = os.path.join(cwd, "LOG")

        args.version = __version__
        args.fastqc_ver = fastqc_ver
        args.trim_galore_ver = trim_galore_ver
        args.cutadapt_ver = cutadapt_ver
        args.rseqc_ver = rseqc_ver
        args.bowtie2_ver = bowtie2_ver
        args.picardtools_ver = picardtools_ver
        args.tophat2_ver = tophat2_ver
        args.feature_counts_ver = feature_counts_ver
        args.R_ver = R_ver
        args.samtools_ver = samtools_ver
        args.hisat_ver = hisat_ver
        args.deseq2_ver = deseq2_ver
        args.report_dir = cwd
        if args.paired:
            args.seq_type = "paired-end"
        else:
            args.seq_type = "single-end"

        lines = []
        for arg in dir(args):
            if not arg.startswith("_"):
                if getattr(args, arg):
                    lines.append( arg + "\t" + str(getattr(args, arg)).strip() )
                else:
                    lines.append( arg + "\t" + "NA" )
        with open(os.path.join(cwd,"Report.data"), "w") as f:
            for line in lines:
                line = line.strip()
                if line:
                    f.write(line+"\n")

        ## DESeq2
        if args.sample_info:
            jobs = ["convert -density 200 {pdf} -flatten 1.png".format(pdf=os.path.join(args.main_outdir, "DESeq2", "Fig2.MA_plot.pdf")),
                    "convert -density 200 {pdf} -flatten 2.png".format(pdf=os.path.join(args.main_outdir, "DESeq2", "Fig3.Vulcano_plot.pdf")),
                    "montage {png1} {png2} -geometry +0.0+0.0 -tile 2x1 {output}".format(png1=os.path.join(cwd,"1.png"), png2=os.path.join(cwd,"2.png"), output=os.path.join(cwd,"plots1.png")),
                    "convert -density 200 {pdf} -flatten 3.png".format(pdf=os.path.join(args.main_outdir, "DESeq2", "Fig5.Heatmap.pdf")),
                    "convert -density 200 {pdf} -flatten 4.png".format(pdf=os.path.join(args.main_outdir, "DESeq2", "Fig6.PCA.pdf")),
                    "montage {png1} {png2} -geometry +0.0+0.0 -tile 2x1 {output}".format(png1=os.path.join(cwd,"3.png"), png2=os.path.join(cwd,"4.png"), output=os.path.join(cwd,"plots2.png")),
                    "rm {} {} {} {}".format(os.path.join(cwd, "1.png"),
                               os.path.join(cwd, "2.png"),
                               os.path.join(cwd, "3.png"),
                               os.path.join(cwd, "4.png"),
                               ),
                    ]
            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
            q.join()
            if is_error:
                exit(is_error)

        jobs = ["[ -f {} ] || ( {} cat {}rna-seq-qc/Report_table.R | {}R --vanilla --quiet --args {} {} )".format(os.path.join(cwd,"Report.tsv"), R_libraries_export, script_path, R_path, args.main_indir, args.main_outdir),
                "{} cat {}rna-seq-qc/Rnw2PDF.R | {}R --vanilla --quiet --args {} {}".format(R_libraries_export, script_path, R_path, cwd, os.path.join(script_path,"rna-seq-qc","Report.Rnw")),
                "rm {} {} {} {}".format(os.path.join(cwd, "Report.aux"),
                                        os.path.join(cwd, "Report.tex"),
                                        os.path.join(cwd, "Report.log"),
                                        os.path.join(cwd, "Report.data"),
                                        ),
                ]
        q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
        q.join()
        if is_error:
            exit(is_error)

        # remove files if
        if args.sample_info:
            jobs = ["rm {} {}".format( os.path.join(cwd, "plots1.png"), os.path.join(cwd, "plots2.png"), ),
                    ]
            q.put(Qjob(jobs, cwd=cwd, logfile=logfile, shell=True, backcopy=True, keep_temp=False))
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
        print "Temp dir:", temp_dir
        print "Genome:", args.genome
        print "Input dir:", args.indir
        print "Host:", socket.gethostname()
        print "Threads per process:", args.threads
        print "Files in parallel:", args.parallel
        print "Trim Galore options:", args.trim_galore_opts
        print "TopHat options:", args.tophat_opts
        print "htseq-count options:", args.htseq_count_opts
        print "Seed (random):", args.seed
        print "FASTA index:", args.fasta_index
        print "Genome index (Bowtie2):", args.genome_index
        print "Transcriptome index (TopHat2):", args.transcriptome_index
        print "Mapping program:", args.mapping_prg
        print "HISAT2 index:", args.hisat_index
        print "Report secondary alignments:", args.report_secondary_alignments
        print "Count program:", args.count_prg
        print "GTF:", args.gtf
        print "BED:", args.bed
        print "Gene names:", args.gene_names
        try:
            print "PATH:", os.environ["PATH"]
        except:
            print ""

    ### Output dir
    ##print "Outdir:", args.outdir
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    os.chdir(args.outdir)

    args.analysis_counter = 0   # Counter for analyzes conducted
    indir = args.indir

    check_for_paired_infiles(args, indir, ".fastq.gz")      # sets the args.paired to True if sequences are paired end

    q = Queue()
    for i in range(args.parallel):
        worker = Thread(target=queue_worker, args=(q, args.verbose))
        worker.setDaemon(True)
        worker.start()

    ## FASTQ downsampling
    t1 = datetime.datetime.now()
    if args.fastq_downsample:
        indir = run_fastq_downsampling(args, q, args.indir)
    ### print "Output folder:", indir
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1

    ## Run FastQC
    t1 = datetime.datetime.now()
    run_fastqc(args, q, indir)
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1

    ## Run Trim Galore!
    if args.trim:
        t1 = datetime.datetime.now()
        indir = run_trim_galore(args, q, indir)
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1

    ## Run FastQC on trimmed reads
    if args.trim:
        ## Run FastQC
        t1 = datetime.datetime.now()
        run_fastqc(args, q, indir, analysis_name="FastQC_on_trimmed_reads")
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1

    ## Run library_type
    t1 = datetime.datetime.now()
    run_library_type(args, q, indir)
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1
    try:
        args.library_type = subprocess.check_output("head -n1 {} | cut -f2".format( os.path.join(args.outdir,"library_type", "library_type.txt") ), shell=True).strip()
    except:
        pass

    if args.paired and args.mapping_prg == 'TopHat2':
        ## Run strand_specificity
        t1 = datetime.datetime.now()
        run_distance_metrics(args, q, indir)
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1

    ## Stop here if requested by user (--no-bam)
    if args.no_bam:
        print "\nPipeline finished because of user option '--no-bam'! rna-seq-qc finished (runtime: {})".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), datetime.datetime.now() - start)
        print "Output stored in: {}\n".format(args.outdir)
        exit(0)

    ## Run TopHat
    if args.mapping_prg == 'TopHat2':
        t1 = datetime.datetime.now()
        bam_dir = run_tophat(args, q, indir)
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1
    elif args.mapping_prg == 'HISAT2':
        t1 = datetime.datetime.now()
        bam_dir = run_hisat2(args, q, indir)
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1

    ## Run htseq-count
    t1 = datetime.datetime.now()
    if args.count_prg == "featureCounts":
        count_dir = run_featureCounts(args, q, bam_dir)
    elif args.count_prg == "htseq-count":
        count_dir = run_htseq_count(args, q, bam_dir)
    else:
        print "Error! Unknown feature counting program:", args.count_prg
        exit(1)
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1

    ## RUN bigwig file creation
    if args.bw:
        t1 = datetime.datetime.now()
        run_bigWig(args, q, bam_dir)
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1

    ## Run DESeq2
    if args.sample_info:
        t1 = datetime.datetime.now()
        run_deseq2(args, q, count_dir)
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1

    ## Run project_summary
    t1 = datetime.datetime.now()
    run_project_report(args, q)
    t2 = datetime.datetime.now()
    print "Duration:", t2-t1

    ## Run RSeQC
    if args.rseqc:
        t1 = datetime.datetime.now()
        run_rseqc(args, q, bam_dir)
        t2 = datetime.datetime.now()
        print "Duration:", t2-t1

    return args.outdir


if __name__ == "__main__":
    start = datetime.datetime.now()
    #print "os.environ['PATH']:    ", os.environ['PATH']
    #print "os.system('echo $PATH): ", os.system("echo $PATH")
    #print "Args:", sys.argv
    outdir = main()
    print "\n{} rna-seq-qc finished (runtime: {})".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), datetime.datetime.now() - start)
    if os.path.isdir(outdir):
        print "Output stored in: {}\n".format(outdir)
    if os.path.isfile(os.path.join(outdir,"project_report","Report.pdf")):
        print "Project report:", os.path.join(outdir,"project_report","Report.pdf")
