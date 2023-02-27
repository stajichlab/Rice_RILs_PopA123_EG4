#!/usr/bin/bash -l

#SBATCH -p short -N 1 -n 4 --mem 48gb --out logs/init.log

module load bwa
module load samtools
DIR=genome
GFF=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.gff.gz
GENOME=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
MRNA=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_rna.fna.gz
PROTEIN=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_protein.faa.gz

mkdir -p $DIR
pushd $DIR
for url in $GFF $GENOME $MRNA $PROTEIN
do
 curl -C - -O $url
done
grep ">" GCF_001433935.1_IRGSP-1.0_genomic.fna | grep chromosome | perl -p -e 's/>(\S+) Oryza sativa Japonica Group cultivar Nipponbare chromosome (\d+),.+/$1,$2/' > chrom_nums.csv
popd


GENOMENAME=$DIR/Nipponbare_IRGSP_1.0
if [ ! -f $GENOME.pac ]; then
	bwa index -p $GENOMENAME $DIR/$(basename $GENOME)
fi
