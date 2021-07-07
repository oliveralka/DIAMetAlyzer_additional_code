withMS-DIAL

Consits of the results based on filtering and comparison. 

Filtering was performed by using by the scripts/convert_and_filter_peak_list.py:

The script filters the MS-DIAL peak list to compounds with MS2 references and
convert the retention time to seconds.

See the "test" directory, for an example:
python ../convert_and_filter_peak_list.py -input_dir "./" -output_dir "./test_output/"


The comparison was performed by using the script scripts/comparison_ground_truth_pyprophet_ms-dial.py:

The script compares the filtered MS-DIAL peak list results to the ground truth and pyprohet (e.g. 5 % FDR) results.

See the "test_comparison" directory, for an example:
python ../comparison_ground_truth_pyprophet_ms-dial.py -input_gd "./comp_p251scadd11_052_pyprophet_005_check.tsv" -input_md "./20210616_export_output_lib_all_PestMix1_8Step1Plasma1SWATH20-50.txt" -output_comparison "./test_comparison_output/20210616_export_output_lib_all_PestMix1_8Step1Plasma1SWATH20-50_comparison_test_output.tsv"

The script "vis_comp_MS-DIAL.Rmd" was used for the visualization of the comparison. 


The "param" directory, contains the spectral library and the experiment file used for the MS-DIAL analysis. 

