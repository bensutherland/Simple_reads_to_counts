#!/bin/bash
# Use ref genome gff to create known and unannotated output gtf for each sample
# Note: Requires sorted .bam 

# Global variables
MAPPED_FOLDER="04_mapped"
REFERENCE_GFF="10_reference/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff"

# User variables
NUM_THREADS="16"

# Run Stringtie on each sorted bam file 
ls -1 $MAPPED_FOLDER/*.sorted.bam |
	perl -pe 's/\.hisat2\.sorted\.bam//' |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name)

          # Run stringtie per sample
	  stringtie -p $NUM_THREADS -G $REFERENCE_GFF -o $i".gtf" -l $label $i".hisat2.sorted.bam" 

done

