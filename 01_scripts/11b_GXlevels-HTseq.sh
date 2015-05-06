#!/bin/bash
# use indexed bam files to generate counts for each individual

# need to create an annotation file that says the entire length of each scaffold is in fact a coding region

REFconvertTOOL="/project/lbernatchez/drobo/users/bensuth/chado_test/chado/bin/gmod_fasta2gff3.pl"

#$REFconvertTOOL \
#	--fasta_dir /project/lbernatchez/drobo/users/bensuth/03_RNASeq/SE-reads_assemble-to-counts/05_trinity_output/Trinity.fasta \
#	--gfffilename Trinity.gff3 \
#	--type CDS \
#	--nosequence


#set environment variables
trimmed_FOLDER="03_trimmed"
# COUNT_FOLDER="09_GXlevels"
REFgff3="./Trinity.gff3"

# use htseq-count
ls -1 $trimmed_FOLDER/*fastq.gz.bam | \
    sort -u | 
    while read i
    do
        echo $i
	htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name $i $REFgff3 > $_htseq_counts.txt
    done
