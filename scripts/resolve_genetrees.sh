# !/bin/bash

directory=$1

for file in $directory/*; do
  python3 randomly_resolve.py -i ${file} -o ${file}.resolved
done
