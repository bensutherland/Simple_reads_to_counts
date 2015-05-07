#!/bin/bash
# use indexed bam files to generate counts for each individual

# point to gmod_fasta2gff2   https://github.com/scottcain/chado_test
# set global variables
REFconvertTOOL="/project/lbernatchez/drobo/users/bensuth/chado_test/chado/bin/gmod_fasta2gff3.pl"
REFERENCE=/project/lbernatchez/drobo/users/bensuth/03_RNASeq/SE-reads_assemble-to-counts/05_trinity_output/Trinity.fasta
ASSEMBLED_FOLDER="05_trinity_output"


# generate .gff from the annotation file indicating that each 'scaffold' is in fact a coding region 
$REFconvertTOOL \
	--fasta_dir $REFERENCE \
	--gfffilename $ASSEMBLED_FOLDER/Trinity.gff3 \
	--type CDS \
	--nosequence


# obtain counts
# set global variables
MAPPED_FOLDER="06_mapped"
COUNT_FOLDER="07_gx_levels"
REFgff3="$ASSEMBLED_FOLDER/Trinity.gff3"

# use htseq-count
ls -1 $MAPPED_FOLDER/*fastq.gz.bam | \
    sort -u | 
    while read i
    do
        echo $i
	htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name $i $REFgff3 > "$i"_htseq_counts.txt
    done

mv $MAPPED_FOLDER/*_htseq_counts.txt $COUNT_FOLDER/
