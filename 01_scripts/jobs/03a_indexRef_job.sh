#!/bin/bash
#$ -N trim
#$ -M bensutherland7@gmail.com 
#$ -m beas
#$ -pe smp 2
#$ -l h_vmem=20G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -S /bin/bash

# used to index reference
REFERENCE=/project/lbernatchez/drobo/users/bensuth/00_resources/Ssa_ASM_3.6.fasta.gz

#unzip fasta
gunzip -c /project/lbernatchez/drobo/users/bensuth/00_resources/Ssa_ASM_3.6.fasta.gz > project/lbernatchez/drobo/users/bensuth/00_resources/Ssa_ASM_3.6.fasta

#Index reference
bwa index $REFERENCE

touch finished.indexing.ref
