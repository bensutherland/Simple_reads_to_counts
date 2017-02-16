#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -J bowtie2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=12:00:00
#SBATCH --mem=10000


# build index for reference

# load modules
module load bowtie/2.1.0

# User variables
REFERENCE="/home/bensuth/00_resources/GCA_000233375.4_ICSASG_v2_genomic_unwrapped.fna.gz"
BASE_NAME="/home/bensuth/00_resources/GCA_000233375.4_ICSASG_v2_genomic_unwrapped"

# Index the reference
bowtie2-build $REFERENCE $BASE_NAME
