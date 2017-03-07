#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH -J samstat 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=12:00:00
#SBATCH --mem=10000

# Load module
module load samtools/1.3

# produces samstat .html files for all files in the $MAPPED_FOLDER
# requirement: must have samstat installed and in path
# need to module load samtools


# Global variables
MAPPED_FOLDER="04_mapped"


# run samstat on files 
ls -1 $MAPPED_FOLDER/*.sam |
	sort -u |
	while read i
	do
	  echo $i
      samstat $i
done
