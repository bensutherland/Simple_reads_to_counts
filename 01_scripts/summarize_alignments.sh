#!/bin/bash
# Summarize alignments from a log file

# Global variables
LOG_FOLDER="20_log_files"

# User variables
LOG_FILE="manual_recording_of_alignments_both_combined_2022-10-11.txt"


# Summarize
cat $LOG_FOLDER/$LOG_FILE | grep -E -A 1 '03_trimmed|overall' - | 
	grep -v 'bam_sort_core' - | 
	grep -vE '^--$' | less

