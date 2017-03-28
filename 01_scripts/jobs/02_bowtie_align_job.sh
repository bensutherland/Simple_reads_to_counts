#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH -J bowtie2 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=200:00:00
#SBATCH --mem=50G

time ./01_scripts/02_bowtie2_aln.sh
