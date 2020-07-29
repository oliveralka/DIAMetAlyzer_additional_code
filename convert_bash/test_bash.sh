#!/bin/bash

in_path="/Volumes/elements/OpenSwath_Metabolomics/20181119_Michael_Katrin/M072_SWATH_unzip/01_Plasma"
out_path="/Volumes/elements/OpenSwath_Metabolomics/20181119_Michael_Katrin/M072_SWATH_unzip/01_Plasma/convert"  

for file in $in_path/*.wiff;
  do
    echo $file 
    echo "$out_path/$(basename $file .wiff).mzML"
  done

