#!/bin/bash

for file in *.sites;
do
  name=$(head -n 1 $file | cut -d $'\t' -f 2)
  awk 'NR % 2 == 0 {print}' $file > $name.txt
  grep -vE "(X)" $name.txt | grep -vE "(N)" | grep -vE "(n)" | tr '[:upper:]' '[:lower:]' > $name_filtered.txt
  rm $name.txt
  mv $name_filtered.txt $name.txt
done 
