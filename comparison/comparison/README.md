These script have been used for the comparison of 
the pyprophet results with the ground truth based on
a manual skyline analysis the pesticide mix. 

Here, first a q-value filtered output is generated. 
# generate q-filtered output
echo "Generatre q-filtered output"
python \
../0_scripts/filter_mscore.py \
-in ./251_052_pyprophet.tsv \
-outdir ./1_filtered/

Then the manual skyline analysis is compared on 
the basis of the assay library with the pyprophet 
results filtered with different q-values. This is 
done to check to build a confusion matrix at each
q-value step (e.g. 1% FDR, 3% FDR, 5% FDR, ...).

TP: Compound was detected in Skyline and Pyprophet (within a RT difference of 5 seconds).
FP: Compound was detected in Pyprophet, but not in skyline.
FN: Compound was detected in Skyline but not in Pyprophet.
TN: Compound was not detected using both analysis methods. 

# compare skyline and openswath output
echo "Create comparison"
python \
../0_scripts/generateComparisonDataWithLibraryAsTemplate_3_loop_ms1ms2.py \
-in_assay ../../input_ref/assaylib_20-50_100.tsv \
-in_sky ../../input_ref/20200302_TransitionResults_20-50_100.csv \
-filtered_dir ./1_filtered/ \
-out_dir ./2_compared/

Afterwards a multiconfusion matrix was constructed from the
results of the comparisons at different q-values. 
This was later used to compare the predicted with the actual FDR (based on the ground truth)
and for the estimation of precision and recall. 

# construct confusion matrix
echo "Construct confusion matrix"
python \
../0_scripts/createMultiConfusionMatrixTable.py \
-input_dir ./2_compared/ \
-out ./251_052_confusion_matrix.tsv


The script can be found in the "0_scripts" directory.
"input_ref" contains the assay library (AssayLibraryGenerator) 
and the skyline TransitionsResults
251_052_pyprophet.tsv is an example for the pyprophet output
251_052_confusion_matrix is an example for the multi confusion matrix output of
the "run_comparison.sh". 
