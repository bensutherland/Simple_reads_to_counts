#!/bin/bash
# merge all sorted bam files (indexed by RG defined earlier) into a single bam and index

MergeProg="/project/lbernatchez/drobo/users/bensuth/programs/MergeSamFiles.jar"
MAPPED_FOLDER="06_mapped"
SNP_FOLDER="08_callSNPs"
# NOTE: make sure to use as specific as possible the name for .bam to not capture the expression .bam when want SNP and visa-versa.

# merge .bam files using picard tools
java -Xmx75g -jar $MergeProg $(printf 'INPUT=%s ' $MAPPED_FOLDER/*_dedup.bam) \
	OUTPUT=$SNP_FOLDER/merged.bam ASSUME_SORTED=TRUE \
	MERGE_SEQUENCE_DICTIONARIES=TRUE

# index merged .bam file
samtools index ./08_callSNPs/merged.bam
