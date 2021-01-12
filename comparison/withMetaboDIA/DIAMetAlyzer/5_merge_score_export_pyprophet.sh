#!/bin/bash

echo "merge"
pyprophet merge --template /Volumes/elements/MetaboDIA/Option_1/s45_mt1_resolve_totoc_02_67/20201129_s45_mt1_resolve_totoc_67_02.pqp --out ./merged.osw *.osw

echo "score"
pyprophet score --in ./merged.osw --out ./scored_ms1ms2.osw --level ms1ms2 --ss_main_score var_isotope_correlation_score --ss_score_filter metabolomics

echo "export"
pyprophet export-compound --in ./scored_ms1ms2.osw --out export_pyprophet_NOFDR.tsv --max_rs_peakgroup_qvalue 1000.0

echo "score_plots"
pyprophet export-compound --in ./scored_ms1ms2.osw --format score_plots

echo "Done" 

