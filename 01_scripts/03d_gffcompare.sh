#!/bin/bash
# Global variables
MAPPED_FOLDER="04_mapped"
REF_FOLDER="10_reference"

# User variables
REFERENCE_GFF="GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff"
NUM_THREADS="3"

# Examine how the transcripts compare with the reference annotation
gffcompare -r $REF_FOLDER/$REFERENCE_GFF -G \
	-o merged $MAPPED_FOLDER/stringtie_merged.gtf 

# Move output to storage folder
mv merged.* $MAPPED_FOLDER/gffcompare/

