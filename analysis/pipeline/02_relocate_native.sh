#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 64 -n 1 --mem 128gb --out logs/relocate2_native.%a.log

module load relocate2
module load bwa
module load samtools
module load bowtie2
module load kent-tools
module load bedtools
module load seqtk

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi

SAMPLES=samples.csv
repeat=$(realpath lib/RiceTE.fa)
genome=$(realpath genome/Nipponbare_IRGSP_1.0)
ref_te=$(realpath genome/Nipponbare_IRGSP_1.0.RM.out)
indir=$(realpath input_dir)
origin=$(realpath input)
outdir=$(realpath relocate2_results_raw)
aligner=blat
size=500

start=`date +%s`

outdir=$outdir
IFS=,
tail -n +2 $SAMPLES | sed -n ${N}p | while read STRAIN FILEBASE
do
   mkdir -p $outdir/$STRAIN
   mkdir -p $SCRATCH/$STRAIN
   
   unset IFS
   for a in $(ls $origin/$FILEBASE)
   do
	   if [[ $(echo "$FILEBASE" | grep -c "\-READ") -ne '0' ]]; then 
		   newname=${STRAIN}_$(basename $a | perl -p -e 's/^(\S+)-READ([12])/$2/')
	   else
		   newname=$(basename $a)
	   fi
	   echo "$a --> $newname"
	   ln -s $a $SCRATCH/$STRAIN/$newname
   done
#   ls -l $SCRATCH/$STRAIN
   relocaTE2.py --te_fasta $repeat --genome_fasta $genome --fq_dir $SCRATCH/$STRAIN --mate_1_id _1 --mate_2_id _2 \
	--outdir $outdir/$STRAIN --reference_ins $ref_te --sample $STRAIN \
	--size $size --step 1234567 --mismatch 2 --cpu $CPU  --aligner $aligner --verbose 4
done
end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

