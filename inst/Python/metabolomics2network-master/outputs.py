# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 17:25:04 2018

@@author: Florence Vinson, Arthur Moreau and Fabien Jourdan

/*******************************************************************************
 * Copyright INRA
 * 
 *  Contact: contact-metexplore@inra.fr
 * 
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *  In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *  The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 ******************************************************************************/
"""
import json

def mapping_to_json(dataset_json_path, metabolites, mapping, output_file, alias_dict):
    """
    Create the output of mapping to json file
    Args:
        arg1 (str): path to the json dataset file
        arg2 (dict): metabolite dictionary (provided by parse functions)
        arg3 (dict): mapping dictionary (cf mapping.py)
        arg4 (str): path to the output file
        arg5 (dict): dictionary of aliases in the json file, obtained from parsing conf.txt
    Returns
    
    """
    with open(dataset_json_path, "r") as dataset_file:
        data = json.load(dataset_file)
        with open(output_file, "w") as outfile:
            for metabolite in data:
                mapped = []
                mol_name = metabolite[alias_dict["name_alias"]]
                if mapping[mol_name] != {}:
                    for met_id in mapping[mol_name]:
                        iddb = metabolites[met_id]["iddb"]
                        distance = mapping[mol_name][met_id]["distance"]
                        mapping_type = mapping[mol_name][met_id]["mapping type"]
                        identifiers = mapping[mol_name][met_id]["identifiers"]
                        path = mapping[mol_name][met_id]["path"]
                        datasetname= mapping[mol_name][met_id]["datasetname"]

                        mapped.append({"path" : path, "datasetname" : datasetname, "idsql" : iddb ,"networkid": met_id, "mapping_type" : mapping_type,"distance": str(distance), "identifiers": identifiers})

                    
                metabolite["mapped"] = mapped

            json.dump(data, outfile, indent =4)
    return


def mapping_to_tsv(mapping, output_file):
    """
    Create the output of mapping to tsv file
    Args:
        arg1 (dict): mapping dictionary (cf mapping.py)
        arg2 (str): path to the output file
    Returns
    """
    with open(output_file, "w") as outfile:
        outfile.write(
            "metabolite name" + "\t" + "mapped on id" + "\t" + "mapping types" + "\t" + "distance" + "\t" + "chebi" + "\t" + "inchi" + "\t" + "smiles" + "\t" + "hmdb" + "\t" + "kegg"
            + "\t" + "pubchem" + "\t" + "lipidmaps" + "\t" + "swisslipids")
        for metabolitename in mapping:
            outfile.write("\n")
            ids_mapped = ""
            mapping_types = ""
            distances = ""
            chebi = ""
            inchi = ""
            kegg = ""
            hmdb = ""
            pubchem = ""
            smiles = ""
            lmid = ""
            swissid = ""

            for network_id in mapping[metabolitename]:
                ids_mapped += network_id + ";"
                distances += str(mapping[metabolitename][network_id]["distance"]) + ";"
                mapping_types += str(mapping[metabolitename][network_id]["mapping type"]) + ";"
                if "chebi" in mapping[metabolitename][network_id]["identifiers"]:
                    chebi += str(mapping[metabolitename][network_id]["identifiers"]["chebi"]) + ";"
                else:
                    chebi += "N/A;"

                if "inchi" in mapping[metabolitename][network_id]["identifiers"]:
                    inchi += str(mapping[metabolitename][network_id]["identifiers"]["inchi"]) + ";"
                else:
                    inchi += "N/A;"

                if "kegg" in mapping[metabolitename][network_id]["identifiers"]:
                    kegg += str(mapping[metabolitename][network_id]["identifiers"]["kegg"]) + ";"
                else:
                    kegg += "N/A;"

                if "hmdb" in mapping[metabolitename][network_id]["identifiers"]:
                    hmdb += str(mapping[metabolitename][network_id]["identifiers"]["hmdb"]) + ";"
                else:
                    hmdb += "N/A;"

                if "pubchem" in mapping[metabolitename][network_id]["identifiers"]:
                    pubchem += str(mapping[metabolitename][network_id]["identifiers"]["pubchem"]) + ";"
                else:
                    pubchem += "N/A;"

                if "smiles" in mapping[metabolitename][network_id]["identifiers"]:
                    smiles += str(mapping[metabolitename][network_id]["identifiers"]["smiles"]) + ";"
                else:
                    smiles += "N/A;"

                if "lipidmaps" in mapping[metabolitename][network_id]["identifiers"]:
                    lmid += str(mapping[metabolitename][network_id]["identifiers"]["lipidmaps"]) + ";"
                else:
                    lmid += "N/A;"

                if "swissid" in mapping[metabolitename][network_id]["identifiers"]:
                    swissid += str(mapping[metabolitename][network_id]["identifiers"]["swissid"]) + ";"
                else:
                    swissid += "N/A;"

            outfile.write(
                metabolitename + "\t" + ids_mapped + "\t" + mapping_types + "\t" + distances + "\t" + chebi + "\t" + inchi
                + "\t" + smiles + "\t" + hmdb + "\t" + kegg + "\t" + pubchem + "\t" + lmid + "\t" + swissid)
    return



