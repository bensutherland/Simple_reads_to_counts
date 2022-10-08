Simple reads to counts

### Disclaimer
This pipeline is made available **with no waranty of usefulness of any kind**.  
Purpose: Multi-map reads to a reference transcriptome or genome to quantified read counts.       

## Overview:
  a) Remove adapters and trim for quality    
  b) Multi-map reads against reference transcriptome or genome    
  c) Use eXpress to estimate gene expression levels per transcript using 'effective counts'  
  
The expression level data can be imported into differential expression analysis software (e.g. edgeR).  

Requires the following:  
`Trimmomatic`   http://www.usadellab.org/cms/?page=trimmomatic  
`bowtie2`       http://bowtie-bio.sourceforge.net/bowtie2/index.shtml        
`samtools`      http://samtools.sourceforge.net    
`eXpress`       https://pachterlab.github.io/eXpress/index.html        
`hisat2`        https://ccb.jhu.edu/software/hisat2/index.shtml    
`stringtie`     https://ccb.jhu.edu/software/stringtie/index.shtml    

## 00. Prepare data

#### Prepare the read data
Put raw fastq.gz single-end data in `02_raw_data`  
Run all scripts from the main directory  

Quality check the data
```
fastqc 02_raw_data/*.fastq.gz -o 02_raw_data/fastqc_raw/ -t 12
multiqc -o 02_raw_data/fastqc_raw/ 02_raw_data/fastqc_raw
```

#### Prepare the assembly
Download the reference genome and copy it to `10_reference`, and prepare for alignment:         
```
# As an example, steps are done with the Pacific oyster genome here
cp -l /scratch2/bsutherland/ref_genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna ./10_reference
cp -l /scratch2/bsutherland/ref_genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff ./10_reference
cp -l /scratch2/bsutherland/ref_genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gtf ./10_reference

# If you have the annotation information, extract splice-site and exon information using HISAT2 package
extract_splice_sites.py ./10_reference/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gtf > 10_reference/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.ss

extract_exons.py ./10_reference/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gtf > 10_reference/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.exon

# Build HISAT2 index. Leave out options for splice-sites or exons if  you do not have a gtf yet
hisat2-build --ss 10_reference/*.ss --exon 10_reference/*.exon 10_reference/GCF_902806645.1_*.fna 10_reference/GCF_902806645.1 

```


### 01. Trim for quality and adaptors
Per sample, trim for quality for SE or PE data. Edit the path to trimmomatic and run:        
Single-end: `01_scripts/01_trimming.sh`     
Paired-end: `01_scripts/01_trimming_PE.sh`      

Quality check the trimmed data     
```
## PE only
# Move singleton output aside, not to be used
mv 03_trimmed/*.single.fq.gz 03_trimmed/singles       

## All
fastqc 03_trimmed/*.paired.fq.gz -o 03_trimmed/fastqc_trimmed -t 12         
multiqc -o 03_trimmed/fastqc_trimmed 03_trimmed/fastqc_trimmed       

```


## Using a reference transcriptome
[if using a reference genome](https://github.com/bensutherland/Simple_reads_to_counts#using-a-reference-genome)       

### 02.A. Multi-map reads against the reference transcriptome     
#### Index
Update REFERENCE variable and index decompressed reference transcriptome with bowtie2.    
`01_scripts/01_bowtie2_build.sh`       
...this only needs to be done once for a reference.      

#### Align
Align each sample from `03_trimmed`, inserting read group IDs.    
Using samtools, convert to .bam, sort, index, and delete .sam and unsorted .bam.    
Single-end: `01_scripts/02_bowtie2_aln_SE.sh`       
Paired-end: `01_scripts/02_bowtie2_aln_PE.sh`        

### 2B) Quantify alignments using eXpress  
Use the sorted bam files to quantify transcript abundances.       
`01_scripts/03_express.sh`      

### 2C) Extract effective counts from eXpress files into edgeR input
Uses files `05_gx_levels/*.xprs`. Open the script `01_scripts/utility_scripts/prepare_gxlevels_matrix.R` in R and use interactively.   
This will output a table entitled `out.matrix.csv`, which can be used as an input to edgeR.    


## Using a reference genome 
[if using a reference transcriptome](https://github.com/bensutherland/Simple_reads_to_counts#using-a-reference-transcriptome)
### 3A) Multi-map reads against a reference genome
First you must unzip the gz assembly to an uncompressed version, as it appears this is required for hisat2-build.    
Once decompressed, index the reference genome with:      
`hisat2-build -p 12 GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna`             

Note: make sure to include the .fna in the 'base name' so that it will remain functional with subsequent scripts.    

Use the script `01_scripts/02_hisat2_aln_PE_to_stringtie.sh` to align paired-end reads to the genome using hisat2, and sort the bam file with samtools. This will require that you update the user variable $REFERENCE with the full path and filename of the reference genome.     
Note: this currently requires that your filenames end in `_R[1|2].paired.fq.gz`        

### 3B) Generate a reference and de novo gff using stringtie 
#### 3B.i) Assemble transcript using a reference genome GTF as a guide
Download the reference genome GTF and GFF from NCBI, place in reference genome folder, and use gunzip to decompress both.     
Update the script `01_scripts/03a_stringtie_bam_to_gtf.sh` to provide the full path to the gff.      
Launch the script, which will use the sorted bam files and create a gtf for each bam file.        
For each sample, stringtie will assemble transcripts based on the reference genome annotation and will find novel transcripts (unannotated) in the reference genome.    

#### 3B.ii) Prepare necessary files for stringtie merge 
Create a `mergelist.txt` and save to `00_archive`. The mergelist will contain each sample's gtf filename (one file per line) using:        
`01_scripts/03b_create_mergelist.sh`       

#### 3B.iii) Merge transcripts from all samples, guided by the reference gff   
Merge sample gtfs into a non-redundant, final merged gtf (i.e., `04_mapped/stringtie_merged.gtf`). Will use the reference genome as a guide to identify novel and known transcripts:        
`01_scripts/03c_stringtie_merge.sh`        

#### 3B.iv) Run gffcompare to generate statistics based on your de novo and reference transcripts
`./01_scripts/03d_gffcompare.sh`    


### 4) Extract read counts from bam files using a gtf
#### 4A) Extract read counts into ctab format with stringtie
Estimate abundances with stringtie:        
`01_scripts/04_stringtie_estimate_abundances.sh`     
Per sample, a folder will be created in `05_gx_levels` with a new .gtf and ctab files.       

#### 4B) Create sample list for extracting expression values
Generate a text file with all sample names and relative paths:      
`01_scripts/05_build_sample_lst.sh`      
...will output to `00_archive/sample_lst.txt`       

#### 4C) Convert gene expression values from ctab format into a count matrix .csv file for edgeR 
Use `prepDE.py` script from stringtie to generate two csv files, which contain the count matrices for genes and transcripts, using coverage values from the output of `stringtie -e` command:      
`prepDE.py -i 00_archive/sample_list.txt -g 05_gx_levels/gene_counts.csv -t 05_gx_levels/transcript_counts.csv --length 150`
note: as needed, update the length flag with the average read length.      



Next: gene expression analysis (edgeR or deseq2 instructions)     
