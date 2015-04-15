#!/bin/bash
#$ -N count
#$ -M bensutherland7@gmail.com
#$ -m beas
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=2:00:00
#$ -cwd
#$ -S /bin/bash

gunzip -c *fastq.gz | grep -c @
