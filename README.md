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
  a) Trim for quality and remove adapters  
  b) Digital normalize libraries to be used for *de novo* reference transcriptome  
  c) Assemble *de novo* reference transcriptome  
  d) Align quality-trimmed sample files against reference  
  e) Obtain expression level raw counts for each contig in reference  
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

# b) Normalization
This step normalizes only those libraries to be used for **de novo** transcriptome assembly  
Edit 01_scripts/02_diginorm.sh by giving path to `insilico_read_normalization.pl` and providing the names of the samples to be used for the assembly, by sample1="your.sample", sample2="your.second.sample" ... etc.  
Samples to be used for normalization will be moved to `04_normalized`, digital normalized, compressed, and concatenated to norm_libs_for_trinity.fq.gz   
requires `insilico_read_normalization.pl`  

Locally:
```
01_scripts/02_diginorm.sh
```

On Katak:
```
qsub 01_scripts/jobs/02_diginorm_job.sh
```

*Tip*: that often the more different individuals included in the assembly, the more contigs are produced. This may be due to alleles dividing the contigs. For this reason, it is suggested to not use all of the individuals for **de novo** assembly, but rather representatives from each condition with the deepest sequencing.

# c) Assemble **de novo** transcriptome
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

This will result in a file called *Trinity.fasta* in your 05_trinity_output folder.

# d) align individual samples against reference

requires `bwa` and `samtools`

### mini-step: index reference with bwa
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

### alignment

Input files are in 03_trimmed/ and output files will be moved to 06_mapped/

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

Endnotes:
You can use the alignment .bam files from step (5) to identify SNPs in your transcripts by following <**new pipeline in development**>
