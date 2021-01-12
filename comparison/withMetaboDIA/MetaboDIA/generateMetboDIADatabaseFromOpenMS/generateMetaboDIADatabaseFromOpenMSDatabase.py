
class ConsensusEntry:
    def __init__(self,
                formula="",
                exact_mass=0.0,
                common_name="",
                systematic_name="",
                db_id="",
                source=""):
        self.formula = formula
        self.exact_mass = exact_mass
        self.common_name = common_name
        self.systematic_name = systematic_name
        self.db_id = db_id
        self.source = source

class StructEntry:
    def __init__(self,
            common_name="",
            systematic_name="",
            db_id=""):
        self.common_name = common_name
        self.systematic_name = systematic_name
        self.db_id = db_id

class MappingEntry:
    def __init__(self,
            formula="",
            exact_mass=0.0,
            db_id="",
            source=""):
        self.formula = formula
        self.exact_mass = exact_mass
        self.db_id = db_id
        self.source = source

class Mapping:
    def __init__(self,
                ids=[],
                entries=[]):
            self.ids = ids
            self.entries = entries

class Struct:
    def __init__(self,
                ids=[],
                entries=[]):
            self.ids = ids
            self.entries = entries

def generateMappingEntries(mapping_paths):
    mapping_entries = {}
    mapping_ids = []
    for path in mapping_paths:   
        f_mapping = open(path,'r')
        for line in f_mapping.readlines():
            one_line = line.strip().split(sep='\t')
            line_count = len(one_line)
            exact_mass = one_line[0]
            formula = one_line[1]
            for i in range(2,line_count):
                entry = MappingEntry()
                entry.formula = formula
                entry.exact_mass = float(exact_mass)
                entry.db_id = one_line[i]
                if entry.db_id.startswith("HMDB"):
                    entry.source = "HMDB"
                else: # try to extract based on the charcters in the database identifier string
                    entry.source = ''.join(filter(lambda i: i.isalpha(), entry.db_id)) 
                mapping_entries[entry.db_id] = entry
                mapping_ids.append(entry.db_id)
        f_mapping.close()
    mapping = Mapping(mapping_ids, mapping_entries)
    return mapping

def generateStructEntries(struct_paths):
    struct_entries = {}
    struct_ids = []
    for path in struct_paths:
        f_struct = open(path,'r')
        for line in f_struct.readlines():
            one_line = line.strip().split(sep='\t')
            entry = StructEntry()
            entry.db_id = one_line[0]
            entry.common_name = one_line[1]
            entry.systematic_name = one_line[1] 
            struct_entries[entry.db_id] = entry
            struct_ids.append(entry.db_id)
        f_struct.close()
    struct = Struct(struct_ids, struct_entries)
    return struct

def generateConsensusEntries(mapping, struct):
    consensus_entries = [] 
    if not list(set(mapping.ids) - set(struct.ids)):
        for element in mapping.ids:
            entry = ConsensusEntry()
            mapping_entry = mapping.entries[element]
            struct_entry = struct.entries[element]
            entry.formula = mapping_entry.formula
            entry.exact_mass = mapping_entry.exact_mass
            entry.common_name = struct_entry.common_name
            entry.systematic_name = struct_entry.systematic_name
            entry.db_id = element
            entry.source = mapping_entry.source
            consensus_entries.append(entry)
    consensus_entries = sorted(consensus_entries, key=lambda entry: entry.exact_mass)
    return consensus_entries

def write(entries, path):
    f = open(path, "w")
    f.writelines("formula" + '\t' + "exactMass" + '\t' + "commonName" + '\t' + "systemicName" + '\t' + "dbID" + '\t' + "source\n")
    for element in entries:
        f.writelines(element.formula + '\t' + str(element.exact_mass) + '\t' + element.common_name + '\t' + element.systematic_name + '\t' + element.db_id + '\t' + element.source + '\n')
    f.close()

if __name__ == "__main__": 

    mapping_paths = ["/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library_metaboDIA/example/HMDBMappingFile.tsv", "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library_metaboDIA/example/lipidmaps_mapping.tsv"]
    struct_paths = ["/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library_metaboDIA/example/HMDB2StructMapping.tsv", "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library_metaboDIA/example/lipismaps_struct.tsv"]
    output_path = "/Volumes/Elements_1/MetaboDIA_data/20200909_data_MetaboDIA/analysis/AMS/library_metaboDIA/example/MetaboDIA_database_oct_2020.tsv"

    mapping = generateMappingEntries(mapping_paths)
    struct = generateStructEntries(struct_paths)
    consensus = generateConsensusEntries(mapping, struct)
    write(consensus, output_path)
    print("--done--")