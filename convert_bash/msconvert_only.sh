#!/bin/bash

inpath=""
outpath=""

for file in $inpath/*.wiff;
  do
    echo $file
    msconvert.exe $file --mzML 
    echo "file - done"
  done
echo "all - DONE" 
