#!/bin/bash
#$ -N bwa
#$ -M $MY_EMAIL_ADDRESS
#$ -m beas
#$ -pe smp 2
#$ -l h_vmem=75G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -S /bin/bash

# point to reference transcriptome
REFERENCE="00_archive/reference_transcriptome/sfontinalis_contigs.fasta.gz"

#unzip fasta
# gunzip -c $REFERENCE.gz > $REFERENCE

#Index reference
bwa index $REFERENCE

touch finished.indexing.ref
