#!/bin/bash

echo "Start: library generation" 

/home/ubuntu/software/openms_build/bin/AssayGeneratorMetabo \
-in /home/ubuntu/data/mzml/*.mzML \
-in_id /home/ubuntu/data/afxml/*.featureXML \
-executable /home/ubuntu/software/sirius/bin/sirius \
-out /home/ubuntu/analysis/s45_mt1_resolve_unknown_totoc_02/20201129_s45_mt1_resolve_unknown_totoc_67_02.tsv \
-fragment_annotation sirius \
-fragment_annotation_score_threshold 0.8 \
-ambiguity_resolution_mz_tolerance 10 \
-ambiguity_resolution_rt_tolerance 20 \
-total_occurrence_filter 0.2 \
-use_known_unknowns \
-decoy_generation \
-decoy_generation_method both \
-use_exact_mass \
-exclude_ms2_precursor \
-min_transitions 3 \
-max_transitions 6 \
-transition_threshold 3.0 \
-preprocessing:filter_by_num_masstraces 1 \
-preprocessing:precursor_mz_tolerance 10 \
-preprocessing:precursor_mz_tolerance_unit ppm \
-preprocessing:precursor_rt_tolerance 5 \
-preprocessing:feature_only \
-sirius:profile qtof \
-sirius:compound_timeout 100 \
-project:processors 28 \
-threads 1
