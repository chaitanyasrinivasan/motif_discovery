#!/bin/bash

file=$1
#name=$(head -n 1 $file | cut -d $'\t' -f 2)
awk 'NR % 2 == 0 {print}' $file > "${file::-4}.txt"
grep -vE "(X)" "${file::-4}.txt" | grep -vE "(N)" | grep -vE "(n)" | tr '[:upper:]' '[:lower:]' > "${file::-4}_filtered.txt"
rm "${file::-4}.txt"
mv "${file::-4}_filtered.txt" "${file::-4}.txt"

