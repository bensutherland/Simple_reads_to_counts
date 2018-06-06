#!/bin/bash
# Obtain gene counts per individual using express

# Set environment variables
MAPPED_FOLDER="04_mapped"
COUNT_FOLDER="05_gx_levels"
EXPRESS="/home/ben/Programs/express-1.5.1-linux_x86_64/express"
REFERENCE="/home/ben/Documents/genomes/GCF_002872995.1_Otsh_v1.0_cds_from_genomic.fna"

# Produce counts per individual with eXpress 
ls -1 $MAPPED_FOLDER/*.sorted.bam |
    sort -u |
    while read i
    do
        echo "Counts for sample" $i
        name=$(basename $i)

        # Run eXpress on the target file
        $EXPRESS $REFERENCE $i

        # Move the result files into the count folder
        mv results.xprs $COUNT_FOLDER/"$name"_results.xprs
        mv params.xprs $COUNT_FOLDER/"$name"_params.xprs 
    done

