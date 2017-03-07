**note: SNP-finding component still in progress**   
**note: some revisions coming for mapping for transcript quantification**    
SE-reads_assemble-to-counts  
Version 0.1  
2015-05-05    

### Disclaimer
This pipeline is made available **with no waranty of usefulness of any kind**.  
Purpose: reference assembly and alignment of single-end data
It is primarily built around valuable tools developed by other groups (see 'Requires' below)  
and components of the IBIS-Trinity-Pipeline: https://github.com/enormandeau/trinity_pipeline_ibis  
and the very useful Simple Fools Guide from the Palumbi lab: http://sfg.stanford.edu/guide.html

# SE-reads_assemble-to-counts
Quality trim single-end data and remove adapters, generate reference transcriptome, map reads to reference transcriptome
## Overview:
  a) Remove adapters and trim for quality    
  b) Multi-map reads against reference transcriptome    
  c) Use Corset to produce clusters for reduction of redundancy of transcriptome and produce counts     
  
The expression level data can be imported into differential expression analysis software (e.g. edgeR).  

Requires the following:  
`Trimmomatic`         http://www.usadellab.org/cms/?page=trimmomatic  
`bowtie2`    
`samtools`            http://samtools.sourceforge.net  
`corset`    

## General comments
Put raw *fastq.gz single-end data in 02_raw_data  
Run all jobs from the main directory  
Job files are specific to Katak at IBIS (slurm), but with some minor editing can be adapted for other servers  

# a) Trim for quality
Generates a fastq file for each library  
requires `Trimmomatic`

Edit 01_scripts/01_trimming.sh by providing the path to `trimmomatic`  
Run locally:
```
01_scripts/01_trimming.sh
```

Run on Katak: 
```
qsub 01_scripts/jobs/01_trimming_job.sh
```

# b) Multi-map reads against the reference transcriptome     

requires `bowtie2` and `samtools`

First index reference with bowtie2 (only need to do once)

Locally:
Ensure reference is in fasta format (not compressed), and change REFERENCE to the path to your reference you want to align against
```
bowtie2-build -f $REFERENCE $REFERENCE
```

On Katak:
```
sbatch 01_scripts/jobs/02_bowtie2_index_job.sh
```

### Alignment

Input files are in 03_trimmed/
Align each sample, inserting read group IDs.    
Using samtools, convert to .bam, sort, index, and remove .sam.    

Locally:
```
01_scripts/02_bowtie2_aln.sh
```

On Katak: 
```
qsub 01_scripts/jobs/02_bowtie_align_job.sh 
```




## Below to be updated (2017-03-07)    

# e) obtain counts for each contig for each individual  
requires `gmod_fasta2gff3.pl` and `htseq-count`

Input files are to be in 06_mapped/

Edit 01_scripts/05_GXlevels.sh by providing the path to gmod_fasta2gff3.pl and the path to the assembled transcriptome to convert it to a .gff3 file for use by `htseq-count`

Locally:
```
01_scripts/05_GXlevels.sh
```

On Katak: 
```
qsub 01_scripts/jobs/05_GXlevels_job.sh
```

The output of the HT-seq script should be ready for input into your preferred analysis pipeline.


# f) SNP finding
Uses code and pipeline from Simple Fools Guide (Palumbi Lab), with useful explanations given for the processes - http://sfg.stanford.edu/SNP.html  
## f) 1) prepare .bam for SNP finding  
Take alignment files generated in step (d) above, and first deduplicate the reads in order to ensure equal coverage for SNP finding. Then merge the .bam files and index.  
requires `MarkDuplicates.jar` from `picard`

Input files are to be in 06_mapped/

Edit 01_scripts/06_dedup-merge-index.sh by providing the path to `MarkDuplicates.jar` and `MergeSamFiles.jar`  

Locally:
```
01_scripts/06_dedup-merge-index.sh
```

On Katak: 
```
qsub 01_scripts/jobs/06_dedup-merge-index_job.sh
```

## f) 2) deal with indels and call high quality SNPs (Training Set)
requires `GenomeAnalysisTK.jar` from `GATK`  

Input files are to be in 08_callSNPs/  
Edit 01_scripts/07_realigner.sh by providing the path to the reference transcriptome and the `GenomeAnalysisTK.jar`   

Locally:
```
01_scripts/07_realigner.sh
```

On Katak: 
```
qsub 01_scripts/jobs/07_realigner_job.sh
```

## f) 3) rediscover SNPs and recalibrate against the training set
requires `GenomeAnalysisTK.jar` from `GATK`

Input files are to be in 08_callSNPs/  
Edit 01_scripts/08_highqualSNPs.sh by providing the path to the reference transcriptome and the `GenomeAnalysisTK.jar`   

Locally:
```
01_scripts/08_highqualSNPs.sh  
```
Then 
```
01_scripts/09_varRecalib.sh  
```



On Katak: 
```
qsub 01_scripts/jobs/08_highqualSNPs_job.sh  
```
On Katak (Then): 
```
qsub 01_scripts/jobs/09_varRecalib_job.sh  
```

Now the final set of high quality SNPs are found in the output file. This can be parsed to obtain genotypes  
