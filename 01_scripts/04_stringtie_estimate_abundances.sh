#!/bin/bash
# Use GTF (created with known and novel features) to quantify expression using -e parameter

# Global variables
MAPPED_FOLDER="04_mapped"
MERGED_GFF="stringtie_merged.gtf" # this is the default name in the pipeline of the created GTF
GXLEVELS_FOLDER="05_gx_levels"

# User variables
NUM_THREADS="12"

# Run StringTie on all samples in mapped folder, identified from generated GTFs 
ls -1 $MAPPED_FOLDER/*.gtf |
        
	# Ignore the merged.gtf in sample name identification
	grep -vE 'stringtie_merged.gtf' |
	perl -pe 's/\.gtf//' |
	sort -u |
	while read i
	do
	  echo $i

          # Generate basename
          name=$(basename $i)
          
          # Make folders per sample
          mkdir $GXLEVELS_FOLDER/$name
          
          # Run stringtie per sample to produce output gtf and ballgown tables
	  stringtie -e -B -p $NUM_THREADS -G $MAPPED_FOLDER/$MERGED_GFF -o $GXLEVELS_FOLDER/$name/$name".gtf" $MAPPED_FOLDER/$name".hisat2.sorted.bam"

          # -e only estimate the abundance of given reference transcripts (requires -G)
          # -B enable output of Ballgown table files (in the same directory as the output GTF) (requires -G, -o recommended)
          # -G <guide_gff>   reference annotation to include in the merging (GTF/GFF3)
          # -o <out_gtf>     output file name for the merged transcripts GTF


done

