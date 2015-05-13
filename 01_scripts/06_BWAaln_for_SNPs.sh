#!/bin/bash
# use bwa mem to map trimmed reads to reference library

###NOTE THIS SCRIPT REQUIRES THAT YOUR REFERENCE IS ALREADY INDEXED (bwa index REFERENCE)

###NOTE this script is very similar to the 04_BWAaln.sh except that it uses the digital normalized files rather than the trimmed only files.

# global variables (note: point to REFERENCE)
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="06_mapped"
REFERENCE="05_trinity_output/sfontinalis_contigs.fasta"

#Create an array to hold the names of all our samples
#Later, we can then cycle through each sample using a simple foor loop
SAMPLES[1]=./$TRIMMED_FOLDER/
SAMPLES[2]=./$TRIMMED_FOLDER/
SAMPLES[3]=./$TRIMMED_FOLDER/

RG[1]='@RG\tID:lib06\tSM:lib06\tPL:Illumina'
RG[2]='@RG\tID:lib07\tSM:lib07\tPL:Illumina'
RG[3]='@RG\tID:lib08\tSM:lib08\tPL:Illumina'

#Map the reads (note add 1 2 3 4 ... up to the number samples given in array above)
for i in 1 2 3
do
    sample=${SAMPLES[${i}]}
    #Map the reads
    bwa mem -t 10 -R ${RG[${i}]} $REFERENCE ${sample} > ${sample}.sam
    samtools view -Sb ${sample}.sam > ${sample}.unsorted.bam  #-S = input sam -b = output bam
    samtools sort ${sample}.unsorted.bam ${sample}
	java -Xmx75g -jar ~/programs/MarkDuplicates.jar INPUT=${sample}.bam OUTPUT=${sample}_dedup.bam METRICS_FILE=${sample}_metricsfile MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=True
	samtools index ${sample}_dedup.bam
done

# clean up space
rm ./$TRIMMED_FOLDER/*.sam ./$TRIMMED_FOLDER/*.unsorted.bam
mv ./$TRIMMED_FOLDER/*_dedup.bam ./$TRIMMED_FOLDER/*_dedup.bam.bai ./$MAPPED_FOLDER/
