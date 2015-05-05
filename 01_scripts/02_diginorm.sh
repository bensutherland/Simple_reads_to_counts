#!/bin/bash
# Normalizing by coverage with 'normalize_by_kmer_coverage.pl'

# Global variables
TRIMMED_FOLDER="03_trimmed"
NORMALIZED_FOLDER="04_normalized"
NORMALIZE_PROGRAM="/prg/trinityrnaseq/trinityrnaseq_r20140717/util/insilico_read_normalization.pl"

# Identify files to be used for assembly
SAMPLES[1]="HI.2494.001.Index_2.lib01_R1_trimmed.fastq.gz"
SAMPLES[2]="HI.2494.001.Index_4.lib02_R1_trimmed.fastq.gz"

# Copy samples of interest to normalized folder
rm -r $NORMALIZED_FOLDER/normalized_reads 2> /dev/null
cp -l $TRIMMED_FOLDER/$SAMPLES $NORMALIZED_FOLDER/

# Digital normalization of samples of interest
ls -1 $NORMALIZED_FOLDER/*trimmed.fastq.gz | \
        perl -pe 's/\.gz//' | \
    sort -u | \
    while read i
    do
        echo 
        echo ">>> Treating $i:"
        echo "-----------------------------------------------------"
        echo "$i"
        echo "-----------------------------------------------------"
        echo 

        # Decompressing reads
        echo "Decompressing reads"
        gunzip -c "$i".gz > "$i"

        $NORMALIZE_PROGRAM \
            --seqType fq \
            --JM 75G \
            --max_cov 50 \
            --single "$i" \
            --output $NORMALIZED_FOLDER \
            --CPU 10

        # Removing decompressed reads
        echo "Removing decompressed reads"
        rm "$i"
    done

# Remove temporary 'ok', 'left', and 'right' links
rm $NORMALIZED_FOLDER/*.fq.ok

# Compress normalized files
echo "Compressing normalized files"
gzip $NORMALIZED_FOLDER/*.fq

# Concatenate all files for Trinity
echo "Concatenating all normalized files for Trinity"
cat $NORMALIZED_FOLDER/*.fq.gz > $NORMALIZED_FOLDER/norm_libs_for_trinity.fq.gz
