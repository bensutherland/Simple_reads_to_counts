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
fastqc 02_raw_data/*.fastq.gz -o 02_raw_data/fastqc_raw/ -t 12
multiqc -o 02_raw_data/fastqc_raw/ 02_raw_data/fastqc_raw
```

## 1) Trim for quality
Generates a fastq file for each library  
requires `Trimmomatic`

Edit 01_scripts/01_trimming.sh by providing the path to `trimmomatic`  
Single-end data: `01_scripts/01_trimming.sh`     
Paired-end data: `01_scripts/01_trimming_PE.sh`   

Quality check the output trimmed data    
If paired-end data, move your singleton trimmed reads to a separate directory; these will not be used:        
`mv 03_trimmed/*.single.fq.gz 03_trimmed/singles`       

Then run fastqc on the paired trimmed data:      
`fastqc 03_trimmed/*.paired.fastq.gz -o 03_trimmed/fastqc_trimmed`         
`multiqc -o 03_trimmed/fastqc_trimmed 03_trimmed/fastqc_trimmed`       

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

## 6) Generate a reference and de novo gff using stringtie 
### a) Assemble transcript using a reference genome GTF as a guide
This will use as an input the reference genome GTF, and alignment .bam files. For each sample, stringtie will assemble transcripts based on the reference genome as well as novel transcripts.    
Launch the script `01_scripts/03a_stringtie_bam_to_gtf.sh`.    
Update the REFERENCE_GFF with the gff for the reference genome assembly.    

### b) Prepare necessary files for stringtie merge 
To do this, you need to generate the file `mergelist.txt`, which contains a name of each sample gtf file, one file on each line.     

Use the `01_scripts/03b_create_mergelist.sh` automated script to create your mergelist.txt in the `00_archive` based on the .gtf files in `04_mapped`.     

### c) Merge transcripts from all samples guided by the reference gff   
Use the `01_scripts/03c_stringtie_merge.sh` to merge sample gtfs using the reference genome as a guide to identify novel and known transcripts. This will produce a non-redundant final, merged gtf.  
The output will be called `04_mapped/stringtie_merged.gtf`    

### d) Run gffcompare to generate statistics based on your de novo and reference transcripts
`./01_scripts/03d_gffcompare.sh`    

## 7) Extract read counts from bam files using the new gff
### a) Use stringtie to extract read counts into ctab format
Per sample, create a folder in `05_gx_levels`, within which a new .gtf and ctab files will be produced.   
`01_scripts/04_stringtie_estimate_abundances.sh`     

### b) Create sample list for extracting expression values
Per sample, produce a text file that contains the sample name, and the relative path to the sample.  This will be used in the next step.     
`05_build_sample_lst.sh`

### c) Convert gene expression values from ctab format into a count matrix .csv file for edgeR 
Use `prepDE.py` (see https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual).    
`prepDE.py -i 00_archive/sample_lst.txt -g 05_gx_levels/gene_counts.csv -t 05_gx_levels/transcript_counts.csv --length 150`


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
