These are a summary of scripts, which can be found in the workflow directory of
the DIAMetAlyzer KNIME workflow. 

1_package_check:
Checks if all necessary packages have been installed for the python instance used in KNIME. 

2_prep_passatutto:
Is used to prepare the fragmentation trees from SIRIUS 4.0.1 to be processed in Passatutto. 

a) rename_and_copy_before_rerooting
Read *.dot files in directory, rename with ScanNr and Name form path and copy them to different location for passatuto processing.

b) library2ms
Parser for SIRIUS 4.0.1 trees to fit passatuto format of a previous SIRIUS version. (Unfortunaltey needed for current passatutto version).

c) extractRerootPCH2allFilterFragments
Tool to extract decoy information from rerooted trees (passatuto) and parse them into the original assay library structure from OpenMS::AssayGeneratorMetabo

3_pyprophet_call:
Python script to run pyprophet merge, score and compound_export function, to
allow for pyprophet processing of the OpenSWATH results.

4_pyprophet_valiation:
Python script to run pyprophet merge, score and compound_export function, to
allow for pyprophet processing of the OpenSWATH results. Here, the merged
and scored .osw files are stored with a different name, that they could be used 
for validation purposes later on. 
 