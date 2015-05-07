#!/bin/bash
# use bwa mem to map trimmed reads to reference library

###NOTE THIS SCRIPT REQUIRES THAT YOUR REFERENCE IS ALREADY INDEXED (bwa index REFERENCE)

# global variables (note: point to REFERENCE)
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="06_mapped"
REFERENCE="05_trinity_output/sfontinalis_contigs.fasta"

#Create an array to hold the names of all our samples
#Later, we can then cycle through each sample using a simple foor loop
SAMPLES[1]=./$TRIMMED_FOLDER/HI.2494.001.Index_12.lib06_R1_trimmed.fastq.gz
SAMPLES[2]=./$TRIMMED_FOLDER/HI.2494.001.Index_13.lib07_R1_trimmed.fastq.gz
SAMPLES[3]=./$TRIMMED_FOLDER/HI.2494.001.Index_14.lib08_R1_trimmed.fastq.gz
SAMPLES[4]=./$TRIMMED_FOLDER/HI.2494.001.Index_2.lib01_R1_trimmed.fastq.gz
SAMPLES[5]=./$TRIMMED_FOLDER/HI.2494.001.Index_4.lib02_R1_trimmed.fastq.gz
SAMPLES[6]=./$TRIMMED_FOLDER/HI.2494.001.Index_5.lib03_R1_trimmed.fastq.gz
SAMPLES[7]=./$TRIMMED_FOLDER/HI.2494.001.Index_6.lib04_R1_trimmed.fastq.gz
SAMPLES[8]=./$TRIMMED_FOLDER/HI.2494.001.Index_7.lib05_R1_trimmed.fastq.gz

RG[1]='@RG\tID:lib06\tSM:lib06\tPL:Illumina'
RG[2]='@RG\tID:lib07\tSM:lib07\tPL:Illumina'
RG[3]='@RG\tID:lib08\tSM:lib08\tPL:Illumina'
RG[4]='@RG\tID:lib01\tSM:lib01\tPL:Illumina'
RG[5]='@RG\tID:lib02\tSM:lib02\tPL:Illumina'
RG[6]='@RG\tID:lib03\tSM:lib03\tPL:Illumina'
RG[7]='@RG\tID:lib04\tSM:lib04\tPL:Illumina'
RG[8]='@RG\tID:lib05\tSM:lib05\tPL:Illumina'

#Map the reads (note add 1 2 3 4 ... up to the number samples given in array above)
for i in 1 2 3 4 5 6 7 8
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
