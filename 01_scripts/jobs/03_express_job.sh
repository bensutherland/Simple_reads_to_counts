#!/bin/bash
#SBATCH --cpus-per-task=3
#SBATCH -J express
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=80:00:00
#SBATCH --mem=30G

time ./01_scripts/03_express.sh
