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
`insilico_read_normalization.pl`  http://trinityrnaseq.sourceforge.net/trinity_insilico_normalization.html  
`Trinity`             http://trinityrnaseq.github.io  
`bwa`                 http://bio-bwa.sourceforge.net  
`samtools`            http://samtools.sourceforge.net  
`gmod_fasta2gff3.pl`  https://github.com/scottcain/chado_test  
`htseq-count`         http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html  

## General comments
Put raw *fastq.gz single-end data in 02_raw_data  
Run all jobs from the main directory  
Job files are specific to Katak at IBIS, but with some minor editing can be adapted for other servers  

# Trim for quality
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

# Normalize library/libraries to be used for reference transcriptome assembly by coverage
Uses libraries put in 03_trimmed
requires `insilico_read_normalization.pl`  

Edit 01_scripts/02_diginorm.sh by giving path to `insilico_read_normalization.pl`

Locally:
```
01_scripts/02_diginorm.sh
```

On Katak:
```
qsub 01_scripts/jobs/02_diginorm_job.sh
```

*Tip*: that often the more different individuals included in the assembly, the more contigs are produced. This may be due to alleles dividing the contigs. For this reason, it is suggested to not use all of the individuals for **de novo** assembly, but rather representatives from each condition with the deepest sequencing.

# Assemble **de novo** transcriptome
Uses normalized libraries in '04_normalized'
requires `Trinity`  

Edit 01_scripts/03_trinity.sh by giving path to `Trinity` assembler

Locally:
```
01_scripts/03_trinity.sh
```

On Katak:
```
qsub 01_scripts/jobs/03_trinity_job.sh
```

# index reference with bwa
Note: only need to do this once  
requires `bwa`  

Locally:
ensure reference is in fasta format (not compressed), and change REFERENCE to the path to your reference you want to align against
```
bwa index REFERENCE
```

On Katak:
```
qsub 01_scripts/jobs/03a_indexRef_job.sh
```

# align individual samples against reference

requires `bwa` and `samtools`

Input files are to be in 06_trimmed_for_mapping/ and output files will be moved to 07_mapped/

Edit 01_scripts/04_BWAaln.sh by giving path and names to each sample and RG indexes (e.g. replace short identifier (bolded for clarity)) for each RG at RG[1]='@RG\tID:**lib208**\tSM:**lib208**\tPL:Illumina'

This step will create an alignment for each sample, and insert read group IDs into each header. Then, using samtools will convert to .bam, sort the .bam, and index it for future work.

Locally:
```
01_scripts/04_BWAaln.sh
```

On Katak: 
```
qsub 01_scripts/jobs/04_BWAaln_job.sh
```

#####
# convert sorted .bam back to .sam (now sorted)
requires `samtools`

Locally:
```
01_scripts/11_GXlevels.sh
```

On Katak: 
```
qsub 01_scripts/jobs/11_GXlevels_job.sh
```

# obtain counts for each contig for each individual

requires `gmod_fasta2gff3.pl` and `htseq-count`

Input files are to be in 06_trimmed_for_mapping/

Edit 01_scripts/11b_GXlevels-HTseq.sh by providing the path to gmod_fasta2gff3.pl and the path to the assembled transcriptome to convert it to a .gff3 file for use by `htseq-count`

Locally:
```
01_scripts/11b_GXlevels-HTseq.sh
```

On Katak: 
```
qsub 01_scripts/jobs/11b_GXlevels-HTseq_job.sh
```
