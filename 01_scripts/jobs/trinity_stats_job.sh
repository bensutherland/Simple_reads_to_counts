#!/bin/bash
#$ -N trinstats
#$ -M $MY_EMAIL_ADDRESS
#$ -m beas
#$ -pe smp 3
#$ -l h_vmem=10G
#$ -l h_rt=4:00:00
#$ -cwd
#$ -S /bin/bash


time /prg/trinityrnaseq/trinityrnaseq_r20140717/util/TrinityStats.pl ./05_trinity_output/Trinity.fasta > ./05_trinity_output/trinity_stats_output.txt
