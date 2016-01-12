infile="$1"
outdir="$2"
deeptools_dir="$3"
samtools="$4/samtools"
threads="$5"

## Defaults (may be overwritten by command line arguments!)
[ -n "$3" ] || deeptools_dir="/package/deeptools-2.0.0/bin/"
[ -n "$4" ] || samtools="/package/samtools-1.2/bin/samtools"
[ -n "$5" ] || threads=4

opts="-p $threads --binSize 25"

[ -f $deeptools_dir/activate ] && source $deeptools_dir/activate


###############################################################################

## Run 
function run {
    echo -e "\n$1\n" 2>&1 | tee -a $outdir/LOG
    eval $1 2>&1 | tee -a $outdir/LOG
}

[ -d $outdir ] || mkdir -p $outdir
cd $outdir


###############################################################################

bname=$(basename ${infile%.bam})
##cname=$(echo $infile | cut -d "/" -f 8)
##echo $bname
##echo $cname

## Both strands
run "[ -f $outdir/${bname}.bothstrands.bam ] || $samtools view -b -f 3 -F 256 -q 20 $infile > $outdir/${bname}.bothstrands.bam"
run "[ -f $outdir/${bname}.bothstrands.bam.bai ] || $samtools index $outdir/${bname}.bothstrands.bam"
run "[ -f $outdir/${bname}.bw ] || $deeptools_dir/bamCoverage -b $outdir/${bname}.bothstrands.bam $opts -o $outdir/${bname}.bw"
##total=$($samtools view -c -F256 ${bname}.bothstrands.bam | awk '{print $1 / 2}')
##echo "Total: $total"

## Forward strand
run "[ -f $outdir/${bname}.fwd1.bam ] || $samtools view -b -f 128 -F 16 $outdir/${bname}.bothstrands.bam > $outdir/${bname}.fwd1.bam"
run "[ -f $outdir/${bname}.fwd2.bam ] || $samtools view -b -f 64  -F 32 $outdir/${bname}.bothstrands.bam > $outdir/${bname}.fwd2.bam"
run "[ -f $outdir/${bname}.fwd.bam ] || $samtools merge -f $outdir/${bname}.fwd.bam $outdir/${bname}.fwd1.bam $outdir/${bname}.fwd2.bam"
run "[ -f $outdir/${bname}.fwd.bam.bai ] || $samtools index $outdir/${bname}.fwd.bam"
##numreads=$($samtools view -c -F256 $outdir/${bname}.fwd.bam | awk '{print $1 / 2}')
##echo "numreads: $numreads"
##scaleFactor=$(echo $numreads $total | awk '{print $1 / $2}')
##echo "scaleFactor: $scaleFactor"
##run "rm $outdir/fwd1.bam $outdir/fwd2.bam"
run "[ -f $outdir/${bname}.fwd.bw ] || $deeptools_dir/bamCoverage -b $outdir/${bname}.fwd.bam $opts -o $outdir/${bname}.fwd.bw"

## Reverse strand
run "[ -f $outdir/${bname}.rev1.bam ] || $samtools view -b -f 144 $outdir/${bname}.bothstrands.bam > $outdir/${bname}.rev1.bam"
run "[ -f $outdir/${bname}.rev2.bam ] || $samtools view -b -f 96  $outdir/${bname}.bothstrands.bam > $outdir/${bname}.rev2.bam"
run "[ -f $outdir/${bname}.rev.bam ] || $samtools merge -f $outdir/${bname}.rev.bam $outdir/${bname}.rev1.bam $outdir/${bname}.rev2.bam"
run "[ -f $outdir/${bname}.rev.bam.bai ] || $samtools index $outdir/${bname}.rev.bam"
##numreads=$($samtools view -c -F256 $outdir/${bname}.rev.bam | awk '{print $1 / 2}')
##echo "numreads: $numreads"
##scaleFactor=$(echo $numreads $total | awk '{print $1 / $2}')
##echo "scaleFactor: $scaleFactor"
##run "rm $outdir/rev1.bam $outdir/rev2.bam"
run "[ -f $outdir/${bname}.rev.bw ] || $deeptools_dir/bamCoverage -b $outdir/${bname}.rev.bam $opts -o $outdir/${bname}.rev.bw"

run "rm $outdir/${bname}.*.bam*"
