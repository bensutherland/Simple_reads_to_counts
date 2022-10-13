#!/bin/bash
# Build a file with sample names and paths
#  will eventually be provided to prepDE.py for moving from StringTie to edgeR

# Global variables
MAPPED_FOLDER="04_mapped"
GXLEVELS_FOLDER="05_gx_levels"

# Remove any previous instances of sample_list
rm 00_archive/sample_list.txt

# Use the GTF filenames in the mapped folder to create a sample list
ls -1 $MAPPED_FOLDER/*.gtf |
        grep -vE 'stringtie_merged.gtf' |
	perl -pe 's/\.gtf//' |
	sort -u |
	while read i
	do

          # Generate basename
          name=$(basename $i)

          echo $name $GXLEVELS_FOLDER/$name/$name".gtf" >> 00_archive/sample_list.txt
          
done

