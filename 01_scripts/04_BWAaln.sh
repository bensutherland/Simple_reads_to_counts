#!/bin/bash
# use bwa mem to map trimmed reads to reference library

#Create an array to hold the names of all our samples
#Later, we can then cycle through each sample using a simple foor loop
SAMPLES[1]=./06_trimmed_for_mapping/HI.2494.001.Index_2.lib01_R1_trimmed.fastq.gz
SAMPLES[2]=./06_trimmed_for_mapping/HI.2494.001.Index_4.lib02_R1_trimmed.fastq.gz
SAMPLES[3]=./06_trimmed_for_mapping/HI.2494.002.Index_7.lib101a_R1_trimmed.fastq.gz
SAMPLES[4]=./06_trimmed_for_mapping/HI.2494.002.Index_12.lib102a_R1_trimmed.fastq.gz

RG[1]='@RG\tID:lib01\tSM:lib01\tPL:Illumina'
RG[2]='@RG\tID:lib02\tSM:lib02\tPL:Illumina'
RG[3]='@RG\tID:lib101a\tSM:lib101a\tPL:Illumina'
RG[4]='@RG\tID:lib102a\tSM:lib102a\tPL:Illumina'

#Create a shell variable to store the location of our reference transcriptome
REFERENCE=/project/lbernatchez/drobo/users/bensuth/03_RNASeq/SE-reads_assemble-to-counts/05_trinity_output/Trinity.fasta

#Index reference transcriptome
bwa index $REFERENCE

#Map the reads
for i in 1 2 3 4
do
    sample=${SAMPLES[${i}]}
    #Map the reads
    bwa mem -t 10 -R ${RG[${i}]} $REFERENCE ${sample} > ${sample}.sam
    samtools view -Sb ${sample}.sam > ${sample}.unsorted.bam  #-S = input sam -b = output bam
    samtools sort ${sample}.unsorted.bam ${sample}
    samtools index ${sample}.bam
    # htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name ${sample}.bam Trinity_all_X.gff3 > ${sample}_htseq_counts.txt
done

mv ./06_trimmed_for_mapping/*.bam ./07_mapped/
mv ./06_trimmed_for_mapping/*.sam ./07_mapped/
