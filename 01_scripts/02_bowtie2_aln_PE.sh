#!/bin/bash
# Map paired-end trimmed reads to reference genome with bowtie2 
# Note: Requires that reference is indexed (see README.md) 
# Note: Requires that input fastq files end in .fq.gz
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="04_mapped"
REFERENCE="/scratch2/bsutherland/ref_txomes/all_unigenes_155.fasta"

# User variables
NUM_THREADS="24"

# Map reads and add RG
ls -1 $TRIMMED_FOLDER/*.paired.fq.gz |
	
	# remove the read 1 or 2 suffix
        sed 's/\_R1\.paired\.fq\.gz//g' |
        sed 's/\_R2\.paired\.fq\.gz//g' |

	# map per sample
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name)
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bowtie2 --end-to-end -k 40 --threads $NUM_THREADS --rg-id $ID -x $REFERENCE -1 $i"_R1.paired.fq.gz" -2 $i"_R2.paired.fq.gz" -S $i.bowtie2.sam
	  samtools view -Sb $i.bowtie2.sam > $i.bowtie2.unsorted.bam
	  samtools sort -n -o $i.bowtie2.sorted.bam $i.bowtie2.unsorted.bam
      rm $i.bowtie2.sam
      rm $i.bowtie2.unsorted.bam
done

# Move output bam files into the mapped folder
mv $TRIMMED_FOLDER/*.bam $MAPPED_FOLDER

