import os
import io
import subprocess
import pandas
import glob

path_to_osw_out = flow_variables['OSW_out_dir']
path_to_assay_lib = flow_variables['decoy_assaylib_out'].replace("file:","")
file = os.path.join(path_to_osw_out,"merged_osw.osw")
path = os.path.join(path_to_osw_out,"*.osw")
inputfiles = glob.glob(path)

commands = ["pyprophet", "merge", "--template", path_to_assay_lib, "--out", file]
for element in inputfiles:
	commands.append(element)

#pyprophet
print("merging")
subprocess.call(commands)

print("scoring")
subprocess.call(["pyprophet", "score", "--in", file, "--level", "ms1ms2", "--ss_main_score", "var_xcorr_shape_weighted"])

print("exporting")
		
# no filtering
subprocess.call(["pyprophet", "export-compound", "--in", file, "--out", file + "_pyprophet_NOFDR_ms1ms2.tsv", "--max_rs_peakgroup_qvalue", "1000.0"])

output_table = pandas.DataFrame({'col1': [1, 2], 'col2': [3, 4]})

