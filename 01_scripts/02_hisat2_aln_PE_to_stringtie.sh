#!/bin/bash
# Map trimmed, paired reads to a reference genome with hisat2 
# Note: Requires that reference is indexed (see README.md) 
# Note: Currently requires that input reads are *R[12].paired.fq.gz
# Note: works with samtools htslib versions
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="04_mapped"
REFERENCE="/scratch2/bsutherland/ref_genomes/GCA_009026015.1_ASM902601v1_genomic.fna"

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
	  hisat2 -x $REFERENCE -k 40 --threads $NUM_THREADS --rg-id $ID -1 $i"_R1.paired.fq.gz" -2 $i"_R2.paired.fq.gz" -S $i.hisat2.sam
	  samtools view -Sb $i.hisat2.sam > $i.hisat2.unsorted.bam
	  samtools sort -o $i.hisat2.sorted.bam $i.hisat2.unsorted.bam

	  # Remove the sam file and unsorted bam to save space
          rm $i.hisat2.sam
	  rm $i.hisat2.unsorted.bam

        done

# Move output sorted bam files into the mapped folder
mv $TRIMMED_FOLDER/*.bam $MAPPED_FOLDER

