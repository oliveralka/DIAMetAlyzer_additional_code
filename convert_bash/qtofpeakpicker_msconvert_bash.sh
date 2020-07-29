#!/bin/bash

in_path=""
out_path=""

for file in $in_path(*.wiff;
 do
   echo "centroiding $file"
   echo "output $(basename $file .wiff).mzML"
   qtofpeakpicker.exe --resolution=20000 --area=1 --threshold=1 --smoothwidth=1.1 --in $file --out "$outpath/$(basename $file .wiff).mzML"
   echo "centroiding - DONE"
   echo "msconvert - START"
   msconvert.exe "$outpath/$(basename $file .wiff).mzML" --mzML --filter --filter "threshold count 150 most-intense" --outfile --out "$outpath/$(basename $file .wiff)_msconvert.mzML"
   echo "msconvert - DONE"
 done
echo "all - DONE" 

