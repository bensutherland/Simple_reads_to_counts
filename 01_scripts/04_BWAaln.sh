#!/bin/bash
# use bwa mem to map trimmed reads to reference library

#Create an array to hold the names of all our samples
#Later, we can then cycle through each sample using a simple foor loop
samples[1]=ORE_wt_rep1
samples[2]=ORE_wt_rep2
samples[3]=ORE_sdE3_rep1
samples[4]=ORE_sdE3_rep2

#Create a shell variable to store the location of our reference genome
reference=/mnt/ebs/trinity_output/Trinity_all_X.fasta

#Make sure we are in the right directory
#Let's store all of our mapping results in /mnt/ebs/rnaseq_mapping2/ to make sure we stay organized
#If this directory already exists, thats ok, but files might get overwritten
cd /mnt/ebs
mkdir rnaseq_mapping2
cd rnaseq_mapping2

#Now we can actually do the mapping and counting
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
    sample=${samples[${i}]}
    #Map the reads
    bwa mem ${reference} /mnt/ebs/trimmed_x/${sample}_1_pe /mnt/ebs/trimmed_x/${sample}_2_pe  > ${sample}.sam
    samtools view -Sb ${sample}.sam > ${sample}.unsorted.bam
    samtools sort ${sample}.unsorted.bam ${sample}
    samtools index ${sample}.bam
    htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name ${sample}.bam Trinity_all_X.gff3 > ${sample}_htseq_counts.txt
done

#set global variables
REFERENCE="sfontinalis_contigs.fasta"
READS="HI*_R1_trimmed.fastq.gz"
# OUTPUT="HI.2494.001.Index_12.lib06_R1_trimmed-out"
HI.2494.001.Index_2.lib01_R1_trimmed.fastq.gz
HI.2494.001.Index_4.lib02_R1_trimmed.fastq.gz
HI.2494.002.Index_12.lib102a_R1_trimmed.fastq.gz
HI.2494.002.Index_7.lib101a_R1_trimmed.fastq.gz


#Indexing reference library:  
gunzip -c ./$REFERENCE.gz > ./$REFERENCE
bwa index $REFERENCE

#align read files
bwa mem -t 10 -R '@RG\tID:lib01\tSM:lib01\tPL:Illumina' $REFERENCE HI.2494.001.Index_2.lib01_R1_trimmed.fastq.gz  > lib01.sam
bwa mem -t 10 -R '@RG\tID:lib02\tSM:lib02\tPL:Illumina' $REFERENCE HI.2494.001.Index_4.lib02_R1_trimmed.fastq.gz  > lib02.sam
bwa mem -t 10 -R '@RG\tID:lib101a\tSM:lib101a\tPL:Illumina' $REFERENCE HI.2494.002.Index_7.lib101a_R1_trimmed.fastq.gz > lib101a.sam
bwa mem -t 10 -R '@RG\tID:lib102a\tSM:lib102a\tPL:Illumina' $REFERENCE HI.2494.002.Index_12.lib102a_R1_trimmed.fastq.gz > lib102a.sam



samtools import $REFERENCE.fai $OUTPUT.sam $OUTPUT.unsorted.bam
samtools sort $OUTPUT.unsorted.bam $OUTPUT
samtools index $OUTPUT.bam

# rm ./sfontinalis_contigs.fasta

