#!/bin/bash
#$ -N diginorm
#$ -M bensutherland7@gmail.com
#$ -m beas
#$ -pe smp 10
#$ -l h_vmem=100G
#$ -l h_rt=72:00:00
#$ -cwd
#$ -S /bin/bash

./01_scripts/02_diginorm.sh
