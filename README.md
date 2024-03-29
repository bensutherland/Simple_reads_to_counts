Simple reads to counts

### Disclaimer
This pipeline is made available **with no waranty of usefulness of any kind**.      
Purpose: Multi-map reads to a reference transcriptome or genome to quantify read counts.       
Run all scripts from the main directory.       

## Overview:
  a) Remove adapters and trim for quality;    
  b) Multi-map reads against reference transcriptome or genome;    
  c) Estimate gene expression levels per gene or transcript for differential expression analysis.     
  
The pipeline currently prepares data for differential expression in edgeR, but the current pipeline does not describe the edgeR workflow, as this will be highly experiment specific.       

**Requires the following:**  
`Trimmomatic`   http://www.usadellab.org/cms/?page=trimmomatic  
`samtools`      http://samtools.sourceforge.net    

Genome-based:     
`hisat2`        https://ccb.jhu.edu/software/hisat2/index.shtml    
`stringtie`     https://ccb.jhu.edu/software/stringtie/index.shtml    

Transcriptome based:      
`bowtie2`       http://bowtie-bio.sourceforge.net/bowtie2/index.shtml        
`eXpress`       https://pachterlab.github.io/eXpress/index.html        


### 00. Prepare data
#### Prepare the read data
Put raw fastq.gz data in `02_raw_data`.        

Raw data quality check:       
```
fastqc 02_raw_data/*.fastq.gz -o 02_raw_data/fastqc_raw/ -t 12
multiqc -o 02_raw_data/fastqc_raw/ 02_raw_data/fastqc_raw
```

#### Genome-based: Prepare the assembly
Download the reference genome, if compressed decompress it, and copy it to `10_reference`, and prepare for alignment:         
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

#### Transcriptome-based: prepare the transcriptome
_to be added_        


### 01. Trim for quality and adaptors
Per sample, trim for quality for SE or PE data. Edit the path to trimmomatic and run:        
Single-end: `01_scripts/01_trimming.sh`     
Paired-end: `01_scripts/01_trimming_PE.sh`      

Trimmed data quality check:     
```
## PE only
# Move singleton output aside, not to be used
mv 03_trimmed/*.single.fq.gz 03_trimmed/singles       

## All
fastqc 03_trimmed/*.paired.fq.gz -o 03_trimmed/fastqc_trimmed -t 12         
multiqc -o 03_trimmed/fastqc_trimmed 03_trimmed/fastqc_trimmed       

```


### 02. Align against a reference transcriptome
[OR if using a reference genome](https://github.com/bensutherland/Simple_reads_to_counts#03-align-against-a-reference-genome)

#### 02.A. Multi-map reads against the reference transcriptome     
**Index**          
Update REFERENCE variable and index decompressed reference transcriptome with bowtie2.    
`01_scripts/01_bowtie2_build.sh`       
...this only needs to be done once for a reference.      

**Align**          
Align each sample from `03_trimmed`, inserting read group IDs.    
Using samtools, convert to .bam, sort, index, and delete .sam and unsorted .bam.    
Single-end: `01_scripts/02_bowtie2_aln_SE.sh`       
Paired-end: `01_scripts/02_bowtie2_aln_PE.sh`        

#### 02.B. Quantify alignments using eXpress  
Use the sorted bam files to quantify transcript abundances.       
`01_scripts/03_express.sh`      

#### 02.C. Extract effective counts from eXpress files into edgeR input
Uses files `05_gx_levels/*.xprs`. Open the script `01_scripts/utility_scripts/prepare_gxlevels_matrix.R` in R and use interactively.   
This will output a table entitled `out.matrix.csv`, which can be used as an input to edgeR.    


### 03. Align against a reference genome 
[OR if using reference transcriptome](https://github.com/bensutherland/Simple_reads_to_counts#02-align-against-a-reference-transcriptome)

#### 3A. Multi-map reads against a reference genome
**Index**       
The genome should already be prepared as described [above](https://github.com/bensutherland/Simple_reads_to_counts#prepare-the-assembly).        

**Align**            
Align PE reads to the genome, convert to BAM, and sort:          
`./01_scripts/02_hisat2_aln_PE_to_stringtie.sh`         
_Note: requires filenames in format of `_R[1|2].paired.fq.gz`_        

(experimental currently) If planning to identify novel transcripts (and known), follow section (3B). If planning to only use the known transcripts from a GFF, then use `01_scripts/04_ST_est_abund_known_only.sh`, pointing to the reference GTF or GFF, then skip ahead to section (4B).        

#### 3B. Generate a reference and de novo gff using StringTie 
**Per sample, assemble transcripts based on BAM and ref GFF**         
Note: the supplied GFF (if exists) should already be in your reference folder (see here)[https://github.com/bensutherland/Simple_reads_to_counts/blob/master/README.md#prepare-the-assembly].         

Update the following script to the guiding GFF (if supplied with genome), and launch:         
`01_scripts/03a_stringtie_bam_to_gtf.sh`         

For each sample, stringtie will assemble transcripts based on the reference genome annotation and will find novel transcripts (unannotated) in the reference genome.    


**Prepare mergelist for StringTie**        
Generate `00_archive/mergelist.txt`, which includes each sample's GTF, one per line:          
`01_scripts/03b_create_mergelist.sh`       


**StringTie merge**
Merge reference GFF and per-sample GTFs into a final, non-redundant GTF (i.e., `04_mapped/stringtie_merged.gtf`).      
Update to point to the reference GFF and launch:         
`01_scripts/03c_stringtie_merge.sh`        


**Evaluate transcript/ gene assemblies**
Run gffcompare to generate statistics based on your de novo and reference transcripts
`./01_scripts/03d_gffcompare.sh`    
...the output will go to `04_mapped/gffcompare/merged.*`          


### 04. Extract feature counts from bam files using the GTF
#### 4A. StringTie extract read counts into ctab format
Estimate abundances with StringTie run with -e parameter:        
`01_scripts/04_stringtie_estimate_abundances.sh`     
Per sample, a folder will be created in `05_gx_levels` with a new .gtf and ctab files.       
note: for more information, see [Using StringTie with DESeq2 and edgeR](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#deseq) 


#### 4B. Create sample and path list for extracting expression values
Note: this step is only needed if you used StringTie to identify novel transcripts, and therefore have a gtf for each sample. If you simply used StringTie to extract known transcript abundances, then skip to prepDE.py step below.      

Generate a text file with all sample names and relative paths:      
`01_scripts/05_build_sample_list.sh`      
...will output to `00_archive/sample_list.txt`       


#### 4C. Convert gene expression values from ctab format into a count matrix .csv file for edgeR 
Use `prepDE.py` script from StringTie to generate count matrices for genes and transcripts based on the coverage values from the above:         
If you created a sample list with per-sample GTF to find novel transcripts:         
`prepDE.py -i 00_archive/sample_list.txt -g 05_gx_levels/gene_counts.csv -t 05_gx_levels/transcript_counts.csv --length 150`
note: as needed, update the length flag with the average read length, which can be found from the multiQC output.      

If you did not search for novel transcripts, and therefore do not have per-sample GTFs:        
`prepDE.py -i 05_gx_levels -g 05_gx_levels/gene_counts.csv -t 05_gx_levels/transcript_counts.csv --length 125`           


### Outputs of pipeline
If using a reference genome, the main data outputs to use will be:     
`05_gx_levels/gene_counts.csv`          
`05_gx_levels/transcript_counts.csv`          

Next: gene expression analysis (edgeR or deseq2 instructions)     
