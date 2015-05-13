#!/bin/bash
# merge all sorted bam files (indexed by RG defined earlier) into a single bam and index
# note that here it is essential to differentiate the .gz.bam for SNPs from those for expr.
samtools merge -rh 00_archive/rg.txt ./08_callSNPs/merged.bam ./07_mapped/*.fq.gz.bam
samtools index ./08_callSNPs/merged.bam
