import sys
sys.path.append(flow_variables['context.workflow.absolute-path'])

import os 
import pandas
pandas.options.mode.chained_assignment = None
import tempfile
import subprocess
import rename_and_copy_before_rerooting_noclick
import library2ms_noclick
import extractRerootPCH2allFilterFragments_noclick

# use current os path separator
sep = os.path.sep
libsep = os.pathsep

# create temporary directory 
tmpdir_renamed = tempfile.TemporaryDirectory()
tmpdir_parsed = tempfile.TemporaryDirectory() 
tmpdir_passatutto = tempfile.TemporaryDirectory() 
absolute_path = flow_variables['context.workflow.absolute-path']

# (1) Rename and copy trees from AssayGeneratorMetabo
sirius_trees = flow_variables['path_sirius_workspace'] 
renamed_trees = tmpdir_renamed.name + sep

rename_and_copy_before_rerooting_noclick.main(sirius_trees,renamed_trees)

# (2) Parse trees to passatuto compatible format
parsed_trees = tmpdir_parsed.name + sep

library2ms_noclick.main(renamed_trees, parsed_trees)

# (3) Call java passatutto from python CLI subprocess 
passatuto_lib = absolute_path + sep + "Passatutto" + sep + "lib" + sep + "passatutto.jar"
math_lib = absolute_path + sep + "Passatutto" + sep + "lib" + sep + "commons-math3-3.4.1.jar"
json_lib = absolute_path + sep + "Passatutto" + sep + "lib" + sep + "json-1.0.jar"
trove_lib = absolute_path + sep + "Passatutto" + sep + "lib" + sep + "trove4j-3.0.3.jar"
all_libs = passatuto_lib + libsep + math_lib + libsep + json_lib + libsep + trove_lib

# passatutto call 

subprocess.call(['java', '-cp', all_libs, 'DecoyDatabaseConstruction', '-target',  parsed_trees.strip('\"'), '-out', tmpdir_passatutto.name.strip('\"'), '-method', 'Reroot'])

# (4) Construct the decoys based on the current assay library from the AssayGeneratoMetabo and the rerooted trees +
if os.name == 'nt':
	openmsassaylib = flow_variables['URI-0'].replace("file:/","")
else:
	openmsassaylib = flow_variables['URI-0'].replace("file:","")
rerootedtrees = tmpdir_passatutto.name.strip('\"') + sep
decoyassaylib = flow_variables['decoy_assaylib_out'].strip('\"')

print(openmsassaylib)
print(rerootedtrees)

# add a CH2 shift on the MS1 precursor level! 
addms1shift = False
min_fragment_mz = flow_variables['min_fragment_mz']
max_fragment_mz = flow_variables['max_fragment_mz']
extractRerootPCH2allFilterFragments_noclick.main(openmsassaylib, rerootedtrees, decoyassaylib, addms1shift, min_fragment_mz, max_fragment_mz)

# check if output is available for further use
if (not os.path.exists(decoyassaylib)):
    print("An error occured - the output file does not exist")
    print("Please check the status of the temporary directories: ")
    print("tmpdir_renamed: " + tmpdir_renamed.name)
    print("tmpdir_parsed: " + tmpdir_parsed.name)
    print("tmpdir_passatutto: " + tmpdir_passatutto.name)
else:    
    print("Target-Decoy assay libary was generated successfully")
    # remove the temporary file
    tmpdir_renamed.cleanup()
    tmpdir_parsed.cleanup()
    tmpdir_passatutto.cleanup()
    flow_variables['decoy_assaylib_out'] = "file:" + decoyassaylib
    # empty dataframe has to be created for node to be run successfully
    output_table = pandas.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
