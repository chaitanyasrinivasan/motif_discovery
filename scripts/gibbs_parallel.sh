#!/bin/bash
input=$1
width=$2
mkdir -p "${input::-4}_splits"
#Shuffle data
perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);' < "${input::-4}_shuffled.txt"
#Partition data
lines_per_file=5
split -d --lines=${lines_per_file} --additional-suffix=.txt "${input::-4}_shuffled.txt" "${input::-4}_splits"/split
rm "${input::-4}_shuffled.txt"
#Create job array of data partitions
ls -d "${input::-4}_splits"/*.txt >> jobs.txt
sbatch --wait submit_jobs.sb ${width}
#Merge output matrices and iterate search 
total_lines=$(wc -l <${input})
#python merge.py jobs.txt "${input::-4}" ${total_lines} ${width}

