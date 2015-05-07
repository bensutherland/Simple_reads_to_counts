#!/bin/bash
#$ -N bwa
#$ -M $MY_EMAIL_ADDRESS
#$ -m beas
#$ -pe smp 10
#$ -l h_vmem=100G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -S /bin/bash

# point to reference transcriptome
REFERENCE=05_trinity_output/sfontinalis_contigs.fasta

#unzip fasta
gunzip -c $REFERENCE.gz > $REFERENCE

#Index reference
bwa index $REFERENCE

touch finished.indexing.ref
