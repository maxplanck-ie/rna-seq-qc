infile="$1"
outdir="$2"
deeptools_dir="$3"
samtools="$4/samtools"
threads="$5"

## Defaults (may be overwritten by command line arguments!)
[ -n "$3" ] || deeptools_dir="/package/deeptools-multiheatmapper/bin/"
[ -n "$4" ] || samtools="/package/samtools-1.2/bin/samtools"
[ -n "$5" ] || threads=8

## bamCoverage options
opts="-v -p $threads --normalizeUsingRPKM --binSize 25 --fragmentLength 1"
 
[ -f $deeptools_dir/activate ] && source $deeptools_dir/activate

########################################################################
## Run 
function run {
    echo -e "\n$1\n" 2>&1 | tee -a $outdir/LOG
    eval $1 2>&1 | tee -a $outdir/LOG
}

## Log
function log {
    echo -e "\n$1\n" 2>&1 | tee -a $outdir/LOG
}
########################################################################

cd $outdir

bname=$(basename ${infile%.bam})
##echo ${bname}

## Both strands
run "[ -f ${bname}.bothstrands.bam ] || $samtools view -b -f 3 -F 256 $infile > ${bname}.bothstrands.bam"
run "[ -f ${bname}.bothstrands.bam.bai ] || $samtools index ${bname}.bothstrands.bam"
run "[ -f ${bname}.bw ] || $deeptools_dir/bamCoverage -b ${bname}.bothstrands.bam $opts -o ${bname}.bw"
total=$($samtools view -c -F256 ${bname}.bothstrands.bam | awk '{print $1 / 2}')
log "Total (${bname}): ${total}"

## Forward strand
run "[ -f ${bname}.fwd1.bam ] || $samtools view -b -f 128 -F 16 ${bname}.bothstrands.bam > ${bname}.fwd1.bam"
run "[ -f ${bname}.fwd2.bam ] || $samtools view -b -f 64  -F 32 ${bname}.bothstrands.bam > ${bname}.fwd2.bam"
run "[ -f ${bname}.fwd.bam ] || $samtools merge -f ${bname}.fwd.bam ${bname}.fwd1.bam ${bname}.fwd2.bam"
run "[ -f ${bname}.fwd.bam.bai ] || $samtools index ${bname}.fwd.bam"
numreads=$($samtools view -c -F256 ${bname}.fwd.bam | awk '{print $1 / 2}')
log "Forward (${bname}): ${numreads}"
scaleFactor=$(echo $numreads $total | awk '{print $1 / $2}')
log "scaleFactor (${bname}): ${scaleFactor}"
run "[ -f ${bname}.fwd.bw ] || $deeptools_dir/bamCoverage -b ${bname}.fwd.bam $opts -o ${bname}.fwd.bw --scaleFactor $scaleFactor" 
    
## Reverse strand
run "[ -f ${bname}.rev1.bam ] || $samtools view -b -f 144 ${bname}.bothstrands.bam > ${bname}.rev1.bam"
run "[ -f ${bname}.rev2.bam ] || $samtools view -b -f 96  ${bname}.bothstrands.bam > ${bname}.rev2.bam"
run "[ -f ${bname}.rev.bam ] || $samtools merge -f ${bname}.rev.bam ${bname}.rev1.bam ${bname}.rev2.bam"
run "[ -f ${bname}.rev.bam.bai ] || $samtools index ${bname}.rev.bam"
numreads=$($samtools view -c -F256 ${bname}.rev.bam | awk '{print $1 / 2}')
log "Reverse (${bname}): ${numreads}"
scaleFactor=$(echo $numreads $total | awk '{print $1 / $2}')
log "scaleFactor (${bname}): ${scaleFactor}"
run "[ -f ${bname}.rev.bw ] || $deeptools_dir/bamCoverage -b ${bname}.rev.bam $opts -o ${bname}.rev.bw --scaleFactor $scaleFactor" 
   
run "rm ${bname}.*.bam*"
