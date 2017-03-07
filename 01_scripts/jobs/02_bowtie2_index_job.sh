#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -J bowtie2-index 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=44:00:00
#SBATCH --mem=10G

# load module
module load bowtie/2.1.0

# Set variables
REFERENCE="/home/bensuth/00_resources/sfontinalis_contigs.fasta"

bowtie2-build -f $REFERENCE $REFERENCE

