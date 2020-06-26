#!/bin/bash

for file in ../data/*.sites;
do
  name=$(head -n 1 $file | cut -d $'\t' -f 2)
  awk 'NR % 2 == 0 {print}' $file > ../data/$name.txt
  grep -vE "(X)" ../data/$name.txt | grep -vE "(N)" | grep -vE "(n)" | tr '[:upper:]' '[:lower:]' > ../data/$name_filtered.txt
  rm ../data/$name.txt
  mv ../data/$name_filtered.txt ../data/$name.txt
done
