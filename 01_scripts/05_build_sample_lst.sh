#!/bin/bash
# Use ref genome gtf to create known and unannotated output gtf for each sample
# Note: Requires sorted .bam 
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
MAPPED_FOLDER="04_mapped"
MERGED_GFF="stringtie_merged.gtf"
GXLEVELS_FOLDER="05_gx_levels"


# User variables
NUM_THREADS="6"

# Run Stringtie 
ls -1 $MAPPED_FOLDER/*.gtf |
        grep -vE 'stringtie_merged.gtf' |
	perl -pe 's/\.gtf//' |
	sort -u |
	while read i
	do
	  #echo $i

          # Generate basename
          name=$(basename $i)

          echo $name $GXLEVELS_FOLDER/$name/$name".gtf" >> 00_archive/sample_lst.txt
          
done

