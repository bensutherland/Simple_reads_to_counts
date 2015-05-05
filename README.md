note: pipeline currently in development

SE-reads_assemble-to-counts
Version 0.1  
2015-05-05

### Disclaimer
This pipeline is made available **with no waranty of usefulness of any kind**.
It has been put together to facilitate reference assembly and alignment sample single-end data
This pipeline uses valuable tools developed by other groups (see 'Requires' below)  
and components of the IBIS-Trinity-Pipeline: https://github.com/enormandeau/trinity_pipeline_ibis  
and the Simple Fools Guide from the Palumbi lab: http://sfg.stanford.edu/guide.html

# SE-reads_assemble-to-counts
Quality trim single-end data and remove adapters, generate reference transcriptome, map reads to reference transcriptome
## Overview:
  1) Trim for quality and remove adapters  
  2) Digital normalize libraries to be used for assembly, use to assemble *de novo* reference  
  3) Align each sample short reads to *de novo* reference  
  4) Obtain expression level raw counts for each contig in reference  
Subsequently, the alignment files can be used in companion pipeline *to be named* for SNP discovery,  
or the expression level data can be imported into differential expression analysis software.  

Requires the following:  
`Trimmomatic`         http://www.usadellab.org/cms/?page=trimmomatic  
`Digital Norm`        http://trinityrnaseq.sourceforge.net/trinity_insilico_normalization.html  
`Trinity`             http://trinityrnaseq.github.io  
`bwa`                 http://bio-bwa.sourceforge.net  
`samtools`            http://samtools.sourceforge.net  
`gmod_fasta2gff3.pl`  https://github.com/scottcain/chado_test  
`htseq-count`         http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html  

## General comments
Raw *fastq.gz single-end data in 02_raw_data; run all jobs from the main directory.
Job files are specific to Katak at IBIS, but with some minor editing can be adapted for other servers.
