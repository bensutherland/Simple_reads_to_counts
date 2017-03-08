#!/bin/bash
# Obtain gene counts per individual using corset 

# Set environment variables
MAPPED_FOLDER="04_mapped"
COUNT_FOLDER="05_gx_levels"
CORSET="/home/bensuth/programs/corset-1.06-linux64/corset"

# Produce counts per individual with corset 
$CORSET \
    -g BC,BC,AC,AC \
    $MAPPED_FOLDER/*.sorted.bam


# Move output
mv clusters.txt $COUNT_FOLDER
mv counts.txt $COUNT_FOLDER
