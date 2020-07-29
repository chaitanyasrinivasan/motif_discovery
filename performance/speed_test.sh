#!/bin/bash

#### Chaitanya Srinivasan ####
# This script measures sequential and parallel motif discovery performance
# over subsets of an input FASTA file

cd ../scripts

TIMEFORMAT=%R #use only real elapsed time
FILE=$1
if [ -z "${FILE}" ];
then
  echo "Error: provide an input file"
  exit 0
fi
SIZE=$(wc -l <${FILE})
for (( i=5; i<=$SIZE; ++i ))
do
	head -n $i $FILE > ../data/data.txt
	(time ./find_motif.sh -i $FILE -w 10 -s) 2>> ../performance/sequential_metrics.txt
	(time ./find_motif.sh -i $FILE -w 10 -p) 2>> ../performance/parallel_metrics.txt
done

# make graph
python plot_metrics.py
exit 1
