#!/bin/bash
# Obtain gene counts per individual using indexed sorted bam files

# NOTE: assumes 'htseq-count' is in your path

# Set environment variables
REF_CONVERT_TOOL="/home/bensuth/programs/chado_test/chado/bin/gmod_fasta2gff3.pl"
REFERENCE="/home/bensuth/03_sfon_txome/SE-reads_assemble-to-counts/00_archive/reference_transcriptome/sfontinalis_contigs.fasta.gz"
MAPPED_FOLDER="06_mapped"
COUNT_FOLDER="07_gx_levels"
REF_gff3="$REFERENCE".gff3

# Produce .gff from transcriptome that indicates each 'scaffold' is a coding sequence
$REF_CONVERT_TOOL \
	--fasta_dir $REFERENCE \
	--gfffilename $REFERENCE.gff3 \
	--type CDS \
	--nosequence

# Produce counts per individual with htseq-count
ls -1 $MAPPED_FOLDER/*.fastq.gz.sorted.bam | 
    sort -u | 
    while read i
    do
        echo "Counts for sample" $i
	    htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name $i $REFgff3 > "$i".htseq.txt
    done

# Clean up mapped folder
mv $MAPPED_FOLDER/*.htseq.txt $COUNT_FOLDER/
