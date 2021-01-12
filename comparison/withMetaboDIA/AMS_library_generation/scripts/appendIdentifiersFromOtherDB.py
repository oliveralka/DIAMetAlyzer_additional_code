
import csv

import math

# https://kodify.net/python/math/round-decimals/
def round_decimals_up(number:float, decimals:int=2):
    """
    Returns a value rounded up to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more")
    elif decimals == 0:
        return math.ceil(number)

    factor = 10 ** decimals
    return math.ceil(number * factor) / factor

class Mapping:
	exact_mass = ""
	formula = ""
	identifiers = []

paths_in = ["/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/scripts/test/HMDBMappingFile.tsv", "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/scripts/test/lipidmaps_mapping.tsv"]
path_out = "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/scripts/test/HMDBLipidMapsMapping.tsv"

mappings = []

for paths in paths_in:
	file = open(paths, 'r')
	lines = file.readlines()

	for line in lines:
		if line.startswith("database"):
			continue

		entry = Mapping()
		splitted = line.strip().split("\t")
		entry.exact_mass = str(float(splitted[0]))
		entry.formula = splitted[1]
		entry.identifiers = splitted[2:]
		mappings.append(entry)

	mappings.sort(key=lambda x: float(x.exact_mass))

unique_mappings = dict()
for entry in mappings:
    key = str(round(float(entry.exact_mass))) + "_" + entry.formula
    if key not in unique_mappings:
        unique_mappings[key] = entry
    else:
        current_mapping = unique_mappings[key]
        current_mapping.identifiers.extend(entry.identifiers)
        if (len(current_mapping.exact_mass) < len(entry.exact_mass)):
            current_mapping.exact_mass = entry.exact_mass 
        unique_mappings[key] = current_mapping

outfile = open(path_out, 'w')
for key,value in unique_mappings.items():
    line = value.exact_mass + "\t" + value.formula + "\t" + '\t'.join(value.identifiers) + '\n'
    outfile.write(line)
outfile.close()




