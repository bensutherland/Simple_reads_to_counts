**note: SNP-finding component still in progress**
**some revisions forthcoming for mapping for transcript quantification**
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
  Follow (e) for expression and (f) for SNPs  
  e) Obtain expression level raw counts for each contig in reference  
  f) merge alignments, search for SNPs with elevated stringency (training set), then in discovery mode  
  
The expression level data can be imported into differential expression analysis software (e.g. edgeR).  

Requires the following:  
`Trimmomatic`         http://www.usadellab.org/cms/?page=trimmomatic  
`insilico_read_normalization.pl`  http://trinityrnaseq.sourceforge.net/trinity_insilico_normalization.html  
`Trinity`             http://trinityrnaseq.github.io  
`bwa`                 http://bio-bwa.sourceforge.net  
`samtools`            http://samtools.sourceforge.net  
`gmod_fasta2gff3.pl`  https://github.com/scottcain/chado_test  
`htseq-count`         http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html  
`picard tools`        http://broadinstitute.github.io/picard/  
`gatk`                https://www.broadinstitute.org/gatk/download/   

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
