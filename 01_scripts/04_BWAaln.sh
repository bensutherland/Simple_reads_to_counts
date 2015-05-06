#!/bin/bash
# use bwa mem to map trimmed reads to reference library

###NOTE THIS SCRIPT REQUIRES THAT YOUR REFERENCE IS ALREADY INDEXED (bwa index REFERENCE)

# global variables (note: point to REFERENCE)
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="06_mapped"
REFERENCE=/project/lbernatchez/drobo/users/bensuth/03_RNASeq/SE-reads_assemble-to-counts/05_trinity_output/Trinity.fasta


#Create an array to hold the names of all our samples
#Later, we can then cycle through each sample using a simple foor loop
SAMPLES[1]=./$TRIMMED_FOLDER/HI.2494.001.Index_2.lib01_R1_trimmed.fastq.gz
SAMPLES[2]=./$TRIMMED_FOLDER/HI.2494.001.Index_4.lib02_R1_trimmed.fastq.gz
SAMPLES[3]=./$TRIMMED_FOLDER/HI.2494.002.Index_7.lib101a_R1_trimmed.fastq.gz
SAMPLES[4]=./$TRIMMED_FOLDER/HI.2494.002.Index_12.lib102a_R1_trimmed.fastq.gz

RG[1]='@RG\tID:lib01\tSM:lib01\tPL:Illumina'
RG[2]='@RG\tID:lib02\tSM:lib02\tPL:Illumina'
RG[3]='@RG\tID:lib101a\tSM:lib101a\tPL:Illumina'
RG[4]='@RG\tID:lib102a\tSM:lib102a\tPL:Illumina'


#Map the reads (note add 1 2 3 4 ... up to the number samples given in array above)
for i in 1 2 3 4
do
    sample=${SAMPLES[${i}]}
    #Map the reads
    bwa mem -t 10 -R ${RG[${i}]} $REFERENCE ${sample} > ${sample}.sam
    samtools view -Sb ${sample}.sam > ${sample}.unsorted.bam  #-S = input sam -b = output bam
    samtools sort ${sample}.unsorted.bam ${sample}
    samtools index ${sample}.bam
done

# clean up space
rm ./$TRIMMED_FOLDER/*.sam ./$TRIMMED_FOLDER/*.unsorted.bam
mv ./$TRIMMED_FOLDER/*.bam ./$TRIMMED_FOLDER/*.bam.bai ./$MAPPED_FOLDER/
