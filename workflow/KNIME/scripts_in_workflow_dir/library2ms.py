# author: Oliver Alka
# date: 25.05.2019
#
# Parser for SIRIUS 4.0.1 trees to fit passatuto format of a 
# previous SIRIUS version.
#  

import csv
import os
import pyopenms
import collections
import re
import click


def formatSirius(path, output):
    counter = 0
    molecularFor = {}
    print(path)
    # specify the output directory!
    if not os.path.isdir(output):
        os.makedirs(output)
    filename = os.path.join(output, os.path.basename(path))

    print(os.path.basename(path))
    current_adduct = re.search(r'\+(\w+)\+', os.path.basename(path)).group(1)
    current_adduct = str(current_adduct)
    # print("current_adduct: ", current_adduct)

    with open(filename, "w+") as f:
        f.write("strict digraph {\n")

        v_count = 1
        with open(path) as fd:
            for line in fd:
                if counter >= 4:
                    # node info
                    string_write = ""

                    if "<BR />" in line:
                        tmp = line.split("<BR />")
                        mass = tmp[2]
                        score = tmp[3].split("</FONT>>")[0]
                        vector = "v" + str(v_count)
                        formular = tmp[0].split(" [")[0]
                        formular = formular.replace("\t", "")
                        # print("formular",formular)

                        # remove current adduct
                        f_radduct = pyopenms.EmpiricalFormula(formular)
                        composition = f_radduct.getElementalComposition()
                        new_comp = {}

                        for k, v in composition.items():
                            k = k.decode('utf-8')
                            new_comp[k] = v

                        # print(print("comp: ", new_comp))
                        new_comp[current_adduct] -= 1  # here "x" should be the adduct
                        new_comp = collections.OrderedDict(sorted(new_comp.items()))
                        f_radduct = ''.join('{}{}'.format(key, value) for key, value in new_comp.items())
                        # print("f_radduct: ", f_radduct)

                        # write to a file
                        # in this line use the formula without the adduct for proper rerooting afterwards
                        string_write = "	" + vector + " [label=\"{}\\n{}\\ncE: [10]\\nScore: {}\"];\n".format(f_radduct,
                                                                                                                 mass,
                                                                                                                 score)
                        molecularFor[formular] = vector
                        # print("updirected molecular info:",string_write)
                        f.write(string_write)
                        v_count += 1
                    # directed path
                    elif " -> " in line:
                        # map formular to node id and give molecule
                        temp = line.split(" -> ")
                        v1 = molecularFor.get(temp[0].replace("\t", ""))
                        # print("temp:",temp[1].split(" "))
                        v2 = molecularFor.get(temp[1].split(" ")[0])
                        labelString = temp[1].split(" ")[1]
                        labelString = labelString.split("[label=<")[1]
                        moleculeName = ""

                        if ("</SUB>" in labelString):

                            moleculeList = labelString.split("</SUB>")

                            for i in range(len(moleculeList) - 1):
                                if "<SUB>" in moleculeList[i]:
                                    atom = moleculeList[i].split("<SUB>")[0]
                                    num = moleculeList[i].split("<SUB>")[1]
                                    moleculeName = moleculeName + atom + num
                            n = len(moleculeList) - 1
                            lastatome = moleculeList[n].split(">];")[0]
                            if lastatome != "":
                                # print("lastatome",lastatome)
                                moleculeName = moleculeName + lastatome

                        else:
                            moleculeName = labelString.split(">")[0]

                        # C3H6O7P -> C3HO2P [label=<H<SUB>5</SUB>O<SUB>5</SUB>>];
                        # v41 -> v32 [label="C6H15N"];
                        # print("moleculeName:",moleculeName)

                        string_write = "	{} -> {} [label=\"{}\"];\n".format(v1, v2, moleculeName)
                        # print("directed info:",string_write)
                        f.write(string_write)
                    else:
                        f.write("\n")

                counter += 1
        f.write("}")


@click.command()
@click.option('--input', '-in', type=click.Path(), help='Input of extracted and renamed trees directory', required=True)
@click.option('--output', '-out', type=click.Path(), help='Output directory for parsed trees - ready for rerooting',
              required=True)
def main(input, output):
    # parse dot file from sirius for passatuto2
    # set directory must include the last /

    for filename in os.listdir(input):
        if filename.endswith(".dot"):
            formatSirius(input + filename, output)
            # print(filename)

    print("Done")
    # formatSirius(directory +"1_C4H4N2O2_M+H+.dot")


if __name__ == "__main__":
    main()
