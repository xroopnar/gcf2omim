#!/bin/env bash

# Check the number of arguments provided ($#)
if [ "$#" -ne 2 ]; then
    # Print the required help message to standard output, describing the bedtools operation
    echo "This script performs a bedtools intersect operation and expects the input to be the path to two bed files. It writes the resulting overlaps to output.txt"
    
    # Exit with a status code of 1, indicating an error or usage issue.
    exit 1
fi

BED_A=$1
BED_B=$2

bedtools intersect -a $BED_A -b $BED_B -wao > $OUTPUT_TSV
