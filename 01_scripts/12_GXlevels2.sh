#!/bin/bash
# uses SFG scripts http://sfg.stanford.edu/expression.html 
# to count expression levels in sorted sam files
# provides a summary statistic file 'summarystat.txt'
# and outputs a file for each input SAMPLENAME_counts.txt

# set environment variables
MAPPED_FOLDER="07_mapped"
#COUNTS_FOLDER="09_GXlevels"
PATH_TO_countxpression="00_archive/countxpression.py"
#PATH_TO_ParseExpression2BigTable_counts="00_archive/ParseExpression2BigTable_counts.py"

$PATH_TO_countxpression 20 20 summarystats.txt $MAPPED_FOLDER/*bam.sam

#Make a file that is a list of all your contig names called ContigNames.txt and the first row called ContigName.
#make contig name file
# ContigNames.txt

#take second column of data (number unique reads mapped to each contig) from counts file from each individual
# $PATH_TO_ParseExpression2BigTable_counts  ContigNames.txt CombinedCounts.txt *counts.txt

