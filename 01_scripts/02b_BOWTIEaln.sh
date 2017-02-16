#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH -J bowtie2 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=44:00:00
#SBATCH --mem=70000

# Map trimmed reads to reference genome with bowtie2 
# Note: Requires that reference is indexed (see README.md) 

# Load modules
module load bowtie/2.1.0

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="04_mapped"
REFERENCE="/home/bensuth/00_resources/GCA_000233375.4_ICSASG_v2_genomic_unwrapped.fna.gz"

# User variables
NUM_THREADS="5"

# map reads and add RG
ls -1 $TRIMMED_FOLDER/*.fastq.gz |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name | perl -pe 's/.*(lib\d+[a-zA-Z]?).*/\1/')
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bowtie2 -q --end-to-end -p $NUM_THREADS -x $REFERENCE -U $i -S $i.bowtie2.sam
      #bwa mem -t $NUM_THREADS -R ${ID} $REFERENCE $i > $i.sam
	  #samtools view -Sb $i.sam > $i.unsorted.bam
	  #samtools sort -o $i.sorted.bam $i.unsorted.bam
	  #samtools index $i.sorted.bam
done

# clean up space
mv ./$TRIMMED_FOLDER/*.sam ./$MAPPED_FOLDER/
