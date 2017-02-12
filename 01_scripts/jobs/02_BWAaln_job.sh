#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH -J bwa 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=44:00:00
#SBATCH --mem=100000

time ./01_scripts/02_BWAaln.sh
