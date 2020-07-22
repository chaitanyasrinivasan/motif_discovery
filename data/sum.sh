#!/bin/bash

for file in *.txt;
do
  sum=$(awk '{sum+=length($0)}END{print sum}' $file)
  avg=$(( sum / 5 ))
  echo "${file} : ${avg}"
done
