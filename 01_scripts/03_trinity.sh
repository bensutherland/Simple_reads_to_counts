#!/bin/bash

#### FIRST SET THE FILES YOU WANT TO USE FOR ASSEMBLY OF REFERENCE ###
# only add the sample name, as below it will be assumed that the file can be found within the $NORMALIZED_FOLDER
sample1="example1"
sample2="example2"
sample3="example3"
# add as many as necessary

### Global variables (modify as needed)
NORMALIZED_FOLDER="04_normalized"
ASSEMBLY_FOLDER="05_trinity_output"
ASSEMBLER="/prg/trinityrnaseq/trinityrnaseq_r20140717/Trinity"

# concatenate samples to be used for Trinity
echo "Concatenating all normalized files for Trinity"
cat $NORMALIZED_FOLDER/$sample1 $NORMALIZED_FOLDER/$sample2 $NORMALIZED_FOLDER/$sample3 > $NORMALIZED_FOLDER/norm_libs_for_trinity.fq.gz

# Assembly with Trinity 
# min_kmer_cov 2 improves rapidity of assembly but is not recommended
echo 
echo "################################################"
echo "# Writing output to: '$ASSEMBLY_FOLDER'"
echo "################################################"
echo 

# clean output directory
rm -r $ASSEMBLY_FOLDER 2> /dev/null

# Launch Trinity
$ASSEMBLER \
    --seqType fq \
    --JM 80G \
    --CPU 10 \
    --single $NORMALIZED_FOLDER/norm_libs_for_trinity.fq.gz \
    --min_contig_length 200 \
    --min_kmer_cov 1 \
    --output $ASSEMBLY_FOLDER

# Cleanup assembly space (optional)
mv $ASSEMBLY_FOLDER/Trinity.fasta .
rm -r $ASSEMBLY_FOLDER/* 2> /dev/null
mv Trinity.fasta $ASSEMBLY_FOLDER
