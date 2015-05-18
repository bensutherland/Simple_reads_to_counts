#!/bin/bash
# uses existing alignments as .bam in mapped folder, marks duplicates, merges and indexes the deduplicated merged .bam

# global variables (note: point to REFERENCE)
MAPPED_FOLDER="06_mapped"
SNP_FOLDER="08_callSNPs"
REFERENCE="05_trinity_output/sfontinalis_contigs.fasta"
markDupProg="/project/lbernatchez/drobo/users/bensuth/programs"
MergeProg="/project/lbernatchez/drobo/users/bensuth/programs/MergeSamFiles.jar"

# dedup, merge, index
ls -1 $MAPPED_FOLDER/*.bam |
	sort -u |
	while read i
	do
	  echo $i
	  java -Xmx75g -jar $markDupProg/MarkDuplicates.jar INPUT="$i" \
		OUTPUT="${i%.bam}"_dedup.bam METRICS_FILE="$i"_metricsfile \
 		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=True
	done

ls -1 $MAPPED_FOLDER/*_dedup.bam |
	sort -u |
	while read i
	do
	  echo $i
	  java -Xmx75g -jar $MergeProg $(printf 'INPUT=%s ' $MAPPED_FOLDER/*_dedup.bam) \
        	OUTPUT=$SNP_FOLDER/merged.bam ASSUME_SORTED=TRUE \
        	MERGE_SEQUENCE_DICTIONARIES=TRUE
	done

samtools index $SNP_FOLDER/merged.bam
