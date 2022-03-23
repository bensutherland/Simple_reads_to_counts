#!/bin/bash
# Index decompressed reference with bowtie2 (only need to do once)

# Global variables
REFERENCE="/scratch2/bsutherland/ref_txomes/all_unigenes_155.fasta"

# User variables
NUM_THREADS="6"

# Index reference
bowtie2-build --threads $NUM_THREADS $REFERENCE $REFERENCE

