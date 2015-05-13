#!/bin/bash
#$ -N HTseq
#$ -M your.email
#$ -m beas
#$ -pe smp 3
#$ -l h_vmem=100G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -S /bin/bash

time ./01_scripts/05_GXlevels.sh

