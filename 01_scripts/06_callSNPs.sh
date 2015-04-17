#!/bin/bash
# first index fasta file, then combine position sorted aln files
# then use mpileup to generate .vcf

#Create a shell variable to store the location of our reference transcriptome
REFERENCE=./05_trinity_output/Trinity.fasta

#index reference transcriptome
samtools faidx $REFERENCE.fai

#combine all RG-added sorted BAM files into one file
cat ./07_mapped/*fastq.gz.bam > ./08_callSNPs/aln.bam

samtools mpileup -uf $REFERENCE.fa ./08_callSNPs/aln.bam | \
	bcftools view -bvcg - > ./08_callSNPs/var.raw.bcf
bcftools view ./08_callSNPs/var.raw.bcf | vcfutils.pl varFilter -D100 > ./08_callSNPs/var.flt.vcf
