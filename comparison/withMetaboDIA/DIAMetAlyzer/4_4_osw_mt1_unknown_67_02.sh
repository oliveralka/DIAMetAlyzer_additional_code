#!/bin/bash

#SBATCH --job-name=67_uk_02
#SBATCH --cpus-per-task 16
#SBATCH --time=48:00:00
#SBATCH --mem=128G

#SBATCH --output=/nfs/wsi/abi/scratch/alka/osw/s45_mt1_resolve_unknown_totoc_67_02/20201205_s45_mt1_resolve_unknown_totoc_67_02.stdout
#SBATCH --error=/nfs/wsi/abi/scratch/alka/osw/s45_mt1_resolve_unknown_totoc_67_02/20201205_s45_mt1_resolve_unknown_totoc_67_02.stderr

echo "Start OSW analysis"

in_path="/nfs/wsi/abi/projects/metabolomics/OpenSWATH/20200909_data_MetaboDIA/analysis/osw/pos_dia"
tmp_cached_mzml_path="/nfs/wsi/abi/scratch/tmp_cached_mzML/s45_mt1_resolve_unknown_totoc_67_02"
library_path="/nfs/wsi/abi/scratch/alka/osw/s45_mt1_resolve_unknown_totoc_67_02/20201129_s45_mt1_resolve_unknown_totoc_67_02.pqp"
out_osw_path="/nfs/wsi/abi/scratch/alka/osw/s45_mt1_resolve_unknown_totoc_67_02/osw"
out_chrom_path="/nfs/wsi/abi/scratch/alka/osw/s45_mt1_resolve_unknown_totoc_67_02/mzml"

for file in /nfs/wsi/abi/projects/metabolomics/OpenSWATH/20200909_data_MetaboDIA/analysis/osw/pos_dia/*.mzML
do

/home/alka/software/openms_branch/openms_branch_s45_resolve/bin/OpenSwathWorkflow \
-in $file \
-tr $library_path \
-out_osw "$out_osw_path/$(basename $file .mzML).osw" \
-out_chrom "$out_chrom_path/$(basename $file .mzML).chrom.mzML" \
-min_upper_edge_dist 1.0 \
-rt_extraction_window 200.0 \
-mz_extraction_window 25.0 \
-mz_extraction_window_ms1 25.0 \
-use_ms1_ion_mobility false \
-readOptions cacheWorkingInMemory \
-tempDirectory $tmp_cached_mzml_path \
-batchSize 250 \
-ms1_isotopes 2 \
-threads 16 \
-Library:retentionTimeInterpretation 'seconds' \
-RTNormalization:lowess:span 0.666666666667 \
-Scoring:rt_normalization_factor 1200.0 \
-Scoring:uis_threshold_sn -1 \
-Scoring:TransitionGroupPicker:use_precursors \
-Scoring:TransitionGroupPicker:compute_peak_shape_metrics \
-Scoring:TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length 7 \
-Scoring:TransitionGroupPicker:PeakIntegrator:integration_type trapezoid \
-Scoring:Scores:use_mi_score false \
-Scoring:Scores:use_ms1_correlation \
-Scoring:Scores:use_ms1_fullscan \
-Scoring:Scores:use_ms1_mi false

done

echo "All done. Have a nice one!" >&2

