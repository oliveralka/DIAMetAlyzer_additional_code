import numpy
import csv

assay_library_header = [
"PrecursorMz",
"ProductMz",
"PrecursorCharge",
"ProductCharge",
"LibraryIntensity",
"NormalizedRetentionTime",
"PeptideSequence",
"ModifiedPeptideSequence",
"PeptideGroupLabel",
"LabelType",
"CompoundName",
"SumFormula",
"SMILES",
"Adducts",
"ProteinId",
"UniprotId",
"GeneName",
"FragmentType",
"FragmentSeriesNumber",
"Annotation",
"CollisionEnergy",
"PrecursorIonMobility",
"TransitionGroupId",
"TransitionId",
"Decoy",
"DetectingTransition",
"IdentifyingTransition",
"QuantifyingTransition",
"Peptidoforms"]

adduct_mass = {"H": 1.007276,
               "Na": 22.989218,
               "K": 38.963158}

class Compound:
    first_formula = ""
    first_adduct = ""
    description = ""
    precursor_mz = ""
    charge = ""
    rt = ""
    fragment_mz = []
    fragment_int = []

class CompoundTransition:
    precursor_mz = 0.0
    product_mz = 0.0
    precursor_charge = "null"
    product_charge = "null"
    library_intensity = 0.0
    normalized_retention_time = "null"
    peptide_sequence = "null"
    modified_peptide_sequence = "null"
    peptide_group_label = "null"
    label_type = "null"
    compound_name = "null"
    sumformula = "null"
    smiles = "null"
    adducts = "null"
    protein_id = "null"
    uniprot_id = "null"	
    gene_name = "null"
    fragment_type = "null"
    fragment_series_number = "null"
    annotation = "null"
    collision_energy = "null"
    precursor_ion_mobility = "null"
    transition_group_id = ""
    transition_id = ""
    decoy = 0
    detecting_transition = 1
    identifying_transition = 0
    quantifying_transition = 1
    peptidoforms = "null"

def read_metabodia_library(in_library_path):
    file = open(in_path, 'r')
    lines = file.readlines()
    file.close()

    # read in compounds
    compoundlist = []
    fragment_mz = []
    fragment_int = []
    for line in lines:
        lastlineoffileemtpy = len(line.strip()) == 0
        if line.startswith("Putative_formula_MS1:") or lastlineoffileemtpy:
            if len(fragment_mz) != 0:
                compound = Compound()
                compound.description = description.strip()
                compound.first_formula = first_formula.strip()
                compound.first_adduct = first_adduct.strip()
                compound.precursor_mz = precursor_mz.strip()
                compound.charge = charge.strip()
                compound.rt = rt.strip()
                compound.fragment_mz = fragment_mz
                compound.fragment_int = fragment_int
                compoundlist.append(compound)
                fragment_mz = []
                fragment_int = []
            if len(line.strip()) != 0:
                description = line.split(" ")[1].strip()
                first_formula = line.split(" ")[1].split("^")[0].strip()
                first_adduct = str("[" + line.split(" ")[1].split("^")[1].split("_")[0].strip() + "]+")
        elif line.startswith("Precursor m/z:"):
            precursor_mz = line.split(" ")[2]
        elif line.startswith("Charge state:"):
            charge = line.split(" ")[2]
        elif line.startswith("RT:"):
            rt = line.split(" ")[1]
        elif line.startswith(tuple('0123456789')):
            fragment_mz.append(line.split("\t")[0].strip())
            fragment_int.append(line.split("\t")[1].strip("\n").strip())
    return compoundlist

in_path = "/Volumes/Samsung_T5/MetaboDIA/script/conversionToAssayLib/test_conversion/20201205_DDA_msconvert_mzXML_10ppm_min_s02_pf00.txt"
out_path = "/Volumes/Samsung_T5/MetaboDIA/script/conversionToAssayLib/test_conversion/DDA_msconvert_mzXML_10ppm_min_s05_pf00.tsv"
min_nr_transitions = 3
max_nr_transitions = 6


filtered_compoundlist = []
compoundlist = read_metabodia_library(in_path)

transitionlist = []
transitions_group_counter = 0

for entry in compoundlist:
    if len(entry.fragment_mz) >= min_nr_transitions:
        filtered_compoundlist.append(entry)

for entry in filtered_compoundlist:
    print(entry.description)

    transition_group_id = entry.description + "_" + entry.rt + "_" + str(transitions_group_counter)

    # sort mass by highest intensity
    mass = numpy.array(entry.fragment_mz)
    intensity = numpy.array(entry.fragment_int)
    ints = intensity.argsort()[::-1][:len(intensity)]
    sort_mass = mass[ints]
    sort_int = intensity[ints] 

    # filter by number of transitions
    entry.fragment_mz = sort_mass[:max_nr_transitions]
    entry.fragment_int = sort_int[:max_nr_transitions]

    # multiple entries for each compound depending on the adducts used! 
    # build compound transitions depending on adduct 
    # one transition for each entry in fragment_mz or fragment_int
    for i in range(0,len(entry.fragment_mz)):
        transition = CompoundTransition()
        transition.precursor_mz =  entry.precursor_mz
        transition.product_mz = entry.fragment_mz[i]
        transition.precursor_charge = entry.charge
        transition.product_charge = "NA"
        transition.library_intensity = entry.fragment_int[i]
        transition.normalized_retention_time = entry.rt
        transition.peptide_sequence = ""
        transition.modified_peptide_sequence = ""
        transition.peptide_group_label = ""
        transition.label_type = ""
        transition.compound_name = entry.description
        transition.sumformula = entry.first_formula
        transition.smiles = "NA"
        transition.adducts = entry.first_adduct
        transition.protein_id = ""
        transition.uniprot_id = ""	
        transition.gene_name = ""
        transition.fragment_type = ""
        transition.fragment_series_number = "-1"
        transition.annotation = "null"
        transition.collision_energy = "-1"
        transition.precursor_ion_mobility = "-1"
        transition.transition_group_id = transition_group_id
        transition.transition_id = transition_group_id + "_" + str(i)
        transition.decoy = 0
        transition.detecting_transition = 1
        transition.identifying_transition = 0
        transition.quantifying_transition = 1
        transition.peptidoforms = ""

        transitionlist.append(transition)

    transitions_group_counter += 1

with open(out_path, 'w', newline='\n') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(assay_library_header)
    for t in transitionlist:
        writer.writerow([
        t.precursor_mz,
        t.product_mz,
        t.precursor_charge,
        t.product_charge,
        t.library_intensity,
        t.normalized_retention_time,
        t.peptide_sequence,
        t.modified_peptide_sequence,
        t.peptide_group_label,
        t.label_type,
        t.compound_name,
        t.sumformula,
        t.smiles,
        t.adducts,
        t.protein_id,
        t.uniprot_id,	
        t.gene_name,
        t.fragment_type,
        t.fragment_series_number,
        t.annotation,
        t.collision_energy,
        t.precursor_ion_mobility,
        t.transition_group_id,
        t.transition_id,
        t.decoy,
        t.detecting_transition,
        t.identifying_transition,
        t.quantifying_transition,
        t.peptidoforms])
