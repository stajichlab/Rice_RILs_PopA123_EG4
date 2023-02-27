#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 64 -n 1 --mem 128gb --out logs/mping_mcclintock_RELOCATE2.%a.log

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

module load mcclintock
SAMPFILE=samples.csv
FASTQFOLDER=input
REF=genome/GCF_001433935.1_IRGSP-1.0_genomic.fna
RM=genome/GCF_001433935.1_IRGSP-1.0_rm.out
RMGFF=genome/GCF_001433935.1_IRGSP-1.0_rm.out.gff
TELIB=lib/mping.fa
OUTDIR=mcclintock_relocate_results
mkdir -p $OUTDIR

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE
do
  LEFT=$(ls $FASTQFOLDER/$FILEBASE | sed -n 1p)
  RIGHT=$(ls $FASTQFOLDER/$FILEBASE | sed -n 2p)
  echo "$LEFT $RIGHT for $FASTQFOLDER/$FILEBASE"
  mcclintock.py -m relocate2 -1 $LEFT -2 $RIGHT -p $CPU \
                -o $OUTDIR/$STRAIN -c $TELIB -r $REF 
done
