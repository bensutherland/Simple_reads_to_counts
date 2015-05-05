#!/bin/bash
#$ -N trim
#$ -M $MY_EMAIL_ADDRESS
#$ -m beas
#$ -pe smp 10
#$ -l h_vmem=100G
#$ -l h_rt=48:00:00
#$ -cwd
#$ -S /bin/bash

./01_scripts/01_trimming.sh
