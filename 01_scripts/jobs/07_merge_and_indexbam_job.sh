#!/bin/bash
#$ -N merge_ind
#$ -M $MY_EMAIL_ADDRESS
#$ -m beas
#$ -pe smp 10
#$ -l h_vmem=100G
#$ -l h_rt=48:00:00
#$ -cwd
#$ -S /bin/bash

time ./01_scripts/07_merge_and_indexbam.sh
