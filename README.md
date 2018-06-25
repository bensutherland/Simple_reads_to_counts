Simple reads to counts

### Disclaimer
This pipeline is made available **with no waranty of usefulness of any kind**.  
Purpose: Multi-mapping alignment of single-end reads to a reference transcriptome and quantification    
It was built within the Bernatchez Lab at IBIS, but is mainly for the authors use and not rigorously tested for others    

## Overview:
  a) Remove adapters and trim for quality    
  b) Multi-map reads against reference transcriptome    
  c) Use eXpress to estimate gene expression levels per transcript using 'effective counts'  
  
The expression level data can be imported into differential expression analysis software (e.g. edgeR).  

Requires the following:  
`Trimmomatic`   http://www.usadellab.org/cms/?page=trimmomatic  
`bowtie2`       http://bowtie-bio.sourceforge.net/bowtie2/index.shtml        
`samtools`      http://samtools.sourceforge.net    
`eXpress`       https://pachterlab.github.io/eXpress/index.html
`hisat2`        https://ccb.jhu.edu/software/hisat2/index.shtml    
`stringtie`     https://ccb.jhu.edu/software/stringtie/index.shtml    

## General comments
Put raw fastq.gz single-end data in `02_raw_data`  
Run all scripts from the main directory  

Quality check the data
```
mkdir 02_raw_data/fastq_raw
fastqc 02_raw_data/*.fq.gz -o 02_raw_data/fastqc_raw/ -t 
multiqc -o 02_raw_data/fastqc_raw/ 02_raw_data/fastqc_raw
```

## 1) Trim for quality
Generates a fastq file for each library  
requires `Trimmomatic`

Edit 01_scripts/01_trimming.sh by providing the path to `trimmomatic`  
Single-end data: `01_scripts/01_trimming.sh`     
Paired-end data: `01_scripts/01_trimming_PE.sh`   

Quality check the output trimmed data    
If paired-end data, make a new directory for the 'single' files, as these will not be used.     
`mkdir 03_trimmed/singles`
`mv 03_trimmed/*.single.* 03_trimmed/singles`

```
mkdir 03_trimmed/fastqc_trimmed
fastqc 03_trimmed/*.paired.fastq.gz -o 03_trimmed/fastqc_trimmed
multiqc -o 03_trimmed/fastqc_trimmed 03_trimmed/fastqc_trimmed
```

## 2) Multi-map reads against the reference transcriptome     
If your data will use a reference genome, skip to step (5).    

Requires `bowtie2` and `samtools`

Index decompressed reference with bowtie2 (only need to do once)
`bowtie2-build --threads 5 -f $REFERENCE $REFERENCE`    


### Alignment

Input files are in `03_trimmed`
Align each sample, inserting read group IDs.    
Using samtools, convert to .bam, sort, index, and remove .sam.    

```
01_scripts/02_bowtie2_aln.sh
```

## 3) Quantify alignments using eXpress  
Uses the sorted bam files to quantify transcript abundances.  
`01_scripts/03_express.sh`  

## 4) Extract effective counts from eXpress files into edgeR input
Uses files `05_gx_levels/*.xprs`. Open the script `01_scripts/utility_scripts/prepare_gxlevels_matrix.R` in R and use interactively.   
This will output a table entitled `out.matrix.csv`, which can be used as an input to edgeR.    

## 5) Multi-map reads against a reference genome
Use the script `02_hisat2_aln_PE_to_stringtie.sh` to align paired-end reads to the genome using hisat2, and with samtools sorting that will be usable for stringtie input.    

## 6) Quantify mappings using stringtie
### Assemble transcript using a reference genome GTF as a guide
Use the script `03_stringtie.sh`. Update the field for the gff file with your gff. 

### Merge transcripts altogether and identify novel and known transcripts
To do this, you need top have `mergelist.txt`, which contains a name of each sample gtf file, one file on each line. 

Use the `01_scripts/03b_create_mergelist.sh` automated script to create your mergelist.txt in the `00_archive`. This will be used in the next step to collect all of your sample gtf files.    

Use the `01_scripts/03c_stringtie_merge.sh` to compare your sample gtf with each other, and with the reference genome gtf to identify novel and known transcripts in a non-redundant final, merged gtf.  





## END 

Untested Section: 
Just in case using alignments against reference genome
**in progress**
Obtain counts for each contig for each individual  
requires `gmod_fasta2gff3.pl` and `htseq-count`

Input files are to be in 06_mapped/

Edit 01_scripts/05_GXlevels.sh by providing the path to gmod_fasta2gff3.pl and the path to the assembled transcriptome to convert it to a .gff3 file for use by `htseq-count`

`01_scripts/05_GXlevels.sh`

The output of the HT-seq script should be ready for input into your preferred analysis pipeline.
