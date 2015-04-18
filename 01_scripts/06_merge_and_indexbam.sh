#!/bin/bash
# merge all sorted bam files using read groups defined earlier into a single bam and index

samtools merge -rh 00_archive/rg.txt ./08_callSNPs/merged.bam ./07_mapped/*.gz.bam
samtools index ./08_callSNPs/merged.bam



