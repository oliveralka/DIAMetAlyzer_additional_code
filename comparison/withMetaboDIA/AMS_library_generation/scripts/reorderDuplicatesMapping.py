
import csv

class Mapping:
    exact_mass = ""
    formula = ""
    identifier = ""
    additional_identifiers = []

path_in = "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/lipidmaps_mapping.tsv"
path_out = "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/lipidmaps_mapping_test_out.tsv"

file = open(path_in, 'r')
lines = file.readlines()

mappings = []

for line in lines:
    if line.startswith("database"):
        continue

    entry = Mapping()
    splitted = line.strip().split("\t")
    entry.exact_mass = splitted[0]
    entry.formula = splitted[1]
    entry.identifier = splitted[2]
    entry.additional_identifiers = []
    mappings.append(entry)

mappings.sort(key=lambda x: x.exact_mass)

unique_mappings = dict()
for entry in mappings:
    key = entry.exact_mass + "_" + entry.formula
    if key not in unique_mappings:
        unique_mappings[key] = entry
    else:
        current_mapping = unique_mappings[key]
        current_mapping.additional_identifiers.append(entry.identifier)
        unique_mappings[key] = current_mapping

outfile = open(path_out, 'w')
for key,value in unique_mappings.items():
    line = value.exact_mass + "\t" + value.formula + "\t" + value.identifier + '\t' + '\t'.join(value.additional_identifiers) + '\n'
    outfile.write(line)
outfile.close()




