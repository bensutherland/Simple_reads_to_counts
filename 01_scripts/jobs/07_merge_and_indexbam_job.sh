#!/bin/bash
#$ -N merge_ind
#$ -M $MY_EMAIL_ADDRESS
#$ -m beas
#$ -pe smp 3
#$ -l h_vmem=100G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -S /bin/bash

time ./01_scripts/07_merge_and_indexbam.sh
