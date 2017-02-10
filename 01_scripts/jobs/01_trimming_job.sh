#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH -J trimming
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=44:00:00
#SBATCH --mem=100000

time ./01_scripts/01_trimming.sh
