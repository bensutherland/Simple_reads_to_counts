#!/bin/bash
# Obtain gene counts per individual using express

# Set environment variables
MAPPED_FOLDER="04_mapped"
COUNT_FOLDER="05_gx_levels"
EXPRESS="/home/bensuth/programs/express-1.5.1-linux_x86_64/express"
REFERENCE="/home/bensuth/00_resources/sfontinalis_contigs.fasta"

# Produce counts per individual with corset 
ls -1 $MAPPED_FOLDER/*.sorted.bam |
    sort -u |
    while read i
    do
        echo "Counts for sample" $i
        name=$(basename $i)
        $EXPRESS $REFERENCE $i
        mv results.xprs $COUNT_FOLDER/"$name"_results.xprs
        mv params.xprs $COUNT_FOLDER/"$name"_params.xprs 
    done

