#!/bin/bash

# currently on branch fix/agm_sirius_index
# to allow for HMDBLipidMaps handling 
# merged to develop + release/2.6.0

files="/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/20200922_preparation_for_ams/single_fxml/*.featureXML"

for file in $files;
do

bn=$(basename $file)
fn="${bn%.*}"

echo "$bn"
echo "$fn"
echo $file

/Users/alka/Documents/work/software/openms_dev/openms_build/bin/AccurateMassSearch \
-in $file  \
-out "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/20200922_AMS/output/mztab/"$fn".mztab" \
-out_annotation "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/20200922_AMS/output/afxml/"$fn".featureXML" \
-db:mapping "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library/HMDBLipidMapsMapping.tsv" \
-db:struct "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library/HMDBLipidMapsStruct.tsv" \
-positive_adducts "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library/restrict_mono_PositiveAdducts.tsv" \
-negative_adducts "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library/restrict_mono_NegativeAdducts.tsv" \
-algorithm:mass_error_value 10 \
-algorithm:mass_error_unit "ppm" \
-algorithm:ionization_mode "positive" \
-algorithm:use_feature_adducts "true" \
-algorithm:mzTab:exportIsotopeIntensities "true"

done


