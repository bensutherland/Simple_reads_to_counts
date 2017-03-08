#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH -J corset
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=80:00:00
#SBATCH --mem=50G

time ./01_scripts/03_corset.sh
