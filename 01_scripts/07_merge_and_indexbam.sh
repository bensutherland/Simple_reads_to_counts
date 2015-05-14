#!/bin/bash
# merge all sorted bam files (indexed by RG defined earlier) into a single bam and index


# NOTE: make sure to use as specific as possible the name for .bam to not capture the expression .bam when want SNP and visa-versa.

samtools merge -rh 00_archive/rg.txt ./08_callSNPs/merged.bam ./06_mapped/*_dedup.bam
samtools index ./08_callSNPs/merged.bam
