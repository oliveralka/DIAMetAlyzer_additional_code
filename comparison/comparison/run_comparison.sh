#!/bin/bash

# only change first "-in" and last "out"

# generate q-filtered output
echo "Generatre q-filtered output"
python \
../0_scripts/filter_mscore.py \
-in ./251_052_pyprophet.tsv \
-outdir ./1_filtered/

# compare skyline and openswath output
echo "Create comparison"
python \
../0_scripts/generateComparisonDataWithLibraryAsTemplate_3_loop_ms1ms2.py \
-in_assay ../../input_ref/assaylib_20-50_100.tsv \
-in_sky ../../input_ref/20200302_TransitionResults_20-50_100.csv \
-filtered_dir ./1_filtered/ \
-out_dir ./2_compared/

# construct confusion matrix
echo "Construct confusion matrix"
python \
../0_scripts/createMultiConfusionMatrixTable.py \
-input_dir ./2_compared/ \
-out ./251_052_confusion_matrix.tsv

echo "Analysis done"
