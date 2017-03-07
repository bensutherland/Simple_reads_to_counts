Simple reads to counts
Version 0.2  
2017-03-07    

### Disclaimer
This pipeline is made available **with no waranty of usefulness of any kind**.  
Purpose: Multi-mapping alignment of single-end reads to a reference transcriptome and quantification    
It was built within the Bernatchez Lab at IBIS, but is mainly for the authors use and not rigorously tested for others    

## Overview:
  a) Remove adapters and trim for quality    
  b) Multi-map reads against reference transcriptome    
  c) Use Corset to produce clusters for reduction of redundancy of transcriptome and produce counts     
  
The expression level data can be imported into differential expression analysis software (e.g. edgeR).  

Requires the following:  
`Trimmomatic`   http://www.usadellab.org/cms/?page=trimmomatic  
`bowtie2`       http://bowtie-bio.sourceforge.net/bowtie2/index.shtml        
`samtools`      http://samtools.sourceforge.net    
`corset`        https://github.com/Oshlack/Corset    

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

