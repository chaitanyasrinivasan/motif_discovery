#!/bin/bash
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
	(time ./gibbs_parallel.sh -i ../data/data.txt -w 10 -t FASTA -s) 2>> ../performance/sequential_metrics.txt
	(time ./gibbs_parallel.sh -i ../data/data.txt -w 10 -t FASTA -p) 2>> ../performance/parallel_metrics.txt
done

# make graph
python plot_metrics.py

exit 1
