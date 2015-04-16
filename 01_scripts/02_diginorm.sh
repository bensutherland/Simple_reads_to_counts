#!/bin/bash
# Normalizing by coverage with 'normalize_by_kmer_coverage.pl'

# Global variables
TRIMMED_FOLDER="03_trimmed"
NORMALIZED_FOLDER="04_normalized"
NORMALIZE_PROGRAM="/prg/trinityrnaseq/trinityrnaseq_r20140717/util/insilico_read_normalization.pl"

# In silico normalization
rm -r $NORMALIZED_FOLDER/normalized_reads 2> /dev/null

ls -1 $TRIMMED_FOLDER/*trimmed.fastq.gz | \
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
