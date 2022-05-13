# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 09:57:01 2018

@author: Florence Vinson, Arthur Moreau and Fabien Jourdan

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


import chebi_graph as chebi_graph
import libchebipy as libchebi



def multimapping(data, network_metabolites):
    """
    Multimapping using every identifier available in the metabolite and lipids dictionaries
    Args:
        arg1 (dict): metabolites, dictionary {met_id : dictionary of identifiers} provided by the parse_metabolites_tsv/json functions in parse module
        arg1 (dict):data, dictionary {lipid_name: dictionary of identifiers} provided by parse_lipidomic_data in parse module
    
    Returns:
        dict: dict: mapping dictionary {lipid: mapped_metabolites:{ distance, mapping type, id_used}}
        
    """
    
    mol_to_metabolites = {}
    
    
    for mol_name in data:
        mol = data[mol_name]
        mapped = {}
        path = []

        for met_id in network_metabolites:
            met = network_metabolites[met_id]
            mapping_identifiers = {}
            

            for inchi in met["inchi"]:
                if inchi != "" and inchi == mol["inchi"]:
                    if "inchi" not in mapping_identifiers:
                        mapping_identifiers["inchi"] = [inchi]
                    else:
                        mapping_identifiers["inchi"].append(inchi)
                                    
            for chebi in met["chebi"]:
                if chebi != "" and (chebi.replace('CHEBI:', '')).isdigit() and chebi == mol["chebi"]:
                    if "chebi" not in mapping_identifiers:
                        mapping_identifiers["chebi"] = [chebi]
                    else:
                        mapping_identifiers["chebi"].append(chebi)
                                        
            for hmdb in met["hmdb"]:
                if hmdb != "" and hmdb == mol["hmdb"]:
                    if "hmdb" not in mapping_identifiers:
                        mapping_identifiers["hmdb"] = [hmdb]
                    else:
                        mapping_identifiers["hmdb"].append(hmdb)
                                    
            for pubchem in met["pubchem"]:
                if pubchem != "" and pubchem == mol["pubchem"]:
                    if "pubchem" not in mapping_identifiers:
                        mapping_identifiers["pubchem"] = [pubchem]
                    else:
                        mapping_identifiers["pubchem"].append(pubchem)
                                        
            for kegg in met["kegg"]:
                if kegg != "" and kegg == mol["kegg"]:
                    if "kegg" not in mapping_identifiers:
                        mapping_identifiers["kegg"] = [kegg]
                    else:
                        mapping_identifiers["kegg"].append(kegg)
                                    
            for smiles in met["smiles"]:
                if smiles != "" and smiles == mol["smiles"] and "[*]" not in smiles: #do not count smiles with * which represent R groups (as we do exact mapping here)
                    if "exact smiles" not in mapping_identifiers:
                        mapping_identifiers["exact smiles"] = [smiles]
                    else:
                        mapping_identifiers["exact smiles"].append(smiles)
                                        
            for lmid in met["lmid"]:
                if lmid != "" and lmid == mol["lmid"]: 
                    if "lipidmaps" not in mapping_identifiers:
                        mapping_identifiers["lipidmaps"] = [lmid]
                    else:
                        mapping_identifiers["lipidmaps"].append(lmid)
                        
            for swissid in met["swissid"]:
                if swissid != "" and swissid == mol["swissid"]:
                    if "swisslipids" not in mapping_identifiers:
                        mapping_identifiers["swisslipids"] = [swissid]
                    else: 
                        mapping_identifiers["swisslipids"].append(swissid)
#                    elif len(mapping_identifiers)["swisslipids"] < 100: #fix a maximum number of id as some of them can match with as much as 6000 swisslipids identifiers, this huge number messes up the lines
#                        mapping_identifiers["swisslipids"].append(swissid)
            
            mapping_details = {"path" : path, "datasetname" : mol_name, "distance" : 0, "mapping type": "exact multimapping","identifiers": mapping_identifiers}
            if mapping_details["identifiers"] != {}:
                mapped[met_id] = mapping_details
        
        mol_to_metabolites[mol_name] = mapped
            
       
    return mol_to_metabolites
 


"""
Functions for chebi class mapping
"""

def convert_to_chebi_entities_once(metabolites):
    """
    fucntion to  convert all metabolites into chebi entities once, because we need them under that from for the mapping
    transofrming them once is much more efficient than doing it at every iteration of the mapping loop in the mapping function
    Args:
        arg1 (dict): dictionary of metabolites provided by parsing functions
    Returns:
        list: list of chebi_entities made with libchebipy
    """

    metabolites_chebi_entities = {}
    
    for met_id in metabolites:
        met = metabolites[met_id]
        metabolites_chebi_entities[met_id] = []
        if met["chebi"] != []:
            for elt in met["chebi"]:
                try:
                    met_chebi_entity = libchebi.ChebiEntity(elt)
                    metabolites_chebi_entities[met_id].append(met_chebi_entity)
                except:
                    print("invalid chebi (" + elt +") for metabolite " + met_id )
    
    return metabolites_chebi_entities
        
    
def chebi_class_mapping(data, metabolites_chebi_entities):
    """
    Mapping using the chebi identifiers and libchebi to buidl a graph and look for parents of the lipids we want to map in the network
    Args:
        arg1 (dict): dictionary of lipids(provided by parsing functions)
        arg2 (list):list of chebi_entities of the metabolites of the network (provided by convert_tochebi_entities_once)
    Returns:
        dict: mapping dictionary {lipid: mapped_metabolites:{ distance, mapping type, id_used}}
    
        
    """
    mol_to_metabolites = {}
    for mol_name in data:
        mapped = {}
        mol = data[mol_name]
        mol_to_metabolites[mol_name] = {}
        mol_chebi = mol["chebi"]
        path = []

        if mol_chebi != "" and (mol_chebi.replace('CHEBI:', '')).isdigit():
            try:
                mol_chebi_entity = libchebi.ChebiEntity(mol_chebi)
                mol_chebi = chebi_graph.get_true_id(mol_chebi_entity)
                # mol_chebi is the chebi in the dataset to be mapped.

                #Create the mapping going upstream
                graph = chebi_graph.create_parent_graph(mol_chebi_entity, {})[1]
                # Create the mapping going downstream
                childen_graph = chebi_graph.create_children_graph(mol_chebi_entity, {})[1]
                mapped = {}
                distance = 99

                #For each metabolites in the network
                for met_id in metabolites_chebi_entities:
                    class_chebi = []
                    for chebi_entity in metabolites_chebi_entities[met_id]:

                        if chebi_entity != None:
                            graph_type= ""
                            met_chebi = chebi_graph.normalize_chebi_id(chebi_graph.get_true_id(chebi_entity))
                            met_conjugates_chebi = chebi_graph.get_conjugate_acid_base(chebi_entity)
                            met_tautomers_chebi = chebi_graph.get_tautomer(chebi_entity)

                            if met_chebi in graph:
                                path = chebi_graph.find_shortest_path(graph, mol_chebi, met_chebi, path=[])
                                new_distance = len(path)-1
                                if new_distance < distance:
                                    distance = new_distance
                                    class_chebi= [met_chebi]
                                    graph_type= "parent"
                                elif new_distance == distance:
                                    class_chebi.append(met_chebi)
                                    graph_type= "parent"

                            if met_chebi in childen_graph:
                                path = chebi_graph.find_shortest_path(childen_graph,  mol_chebi, met_chebi,  path=[])
                                new_distance = len(path) - 1
                                if new_distance < distance:
                                    distance = new_distance
                                    class_chebi = [met_chebi]
                                    graph_type= "child"
                                elif new_distance == distance:
                                    class_chebi.append(met_chebi)
                                    graph_type= "child"

                            for met_conjugate in met_conjugates_chebi:
                                if met_conjugate in graph:
                                    path = chebi_graph.find_shortest_path(graph, mol_chebi, met_conjugate, path=[])
                                    new_distance = len(path)-1 + 0.1
                                    if new_distance < distance:
                                        distance = new_distance
                                        class_chebi= [met_chebi]
                                        graph_type= "parent"
                                    elif new_distance == distance:
                                        class_chebi.append(met_chebi)
                                        graph_type= "parent"

                                if met_conjugate in childen_graph:
                                    path = chebi_graph.find_shortest_path(childen_graph, mol_chebi, met_conjugate, path=[])
                                    new_distance = len(path)-1 + 0.1
                                    if new_distance < distance:
                                        distance = new_distance
                                        class_chebi= [met_chebi]
                                        graph_type= "child"
                                    elif new_distance == distance:
                                        class_chebi.append(met_chebi)
                                        graph_type= "child"

                            for met_tautomer in met_tautomers_chebi:
                                if met_tautomer in graph:
                                    path = chebi_graph.find_shortest_path(graph, mol_chebi, met_tautomer, path=[])
                                    new_distance = len(path)-1 + 0.01
                                    if new_distance < distance:
                                        distance = new_distance
                                        class_chebi= [met_chebi]
                                        graph_type= "parent"
                                    elif new_distance == distance:
                                        class_chebi.append(met_chebi)
                                        graph_type= "parent"

                                if met_tautomer in childen_graph:
                                    path = chebi_graph.find_shortest_path(childen_graph, mol_chebi, met_tautomer, path=[])
                                    new_distance = len(path)-1 + 0.01
                                    if new_distance < distance:
                                        distance = new_distance
                                        class_chebi= [met_chebi]
                                        graph_type= "child"
                                    elif new_distance == distance:
                                        class_chebi.append(met_chebi)
                                        graph_type= "child"

                            result_distance= distance
                            if graph_type == "child":
                                result_distance= -distance

                            mapping_details = {"path" : path, "datasetname" : mol_name, "distance" : result_distance, "mapping type": "chebi class mapping",
                                               "identifiers": {"chebi": class_chebi}}

                            if class_chebi != []:
                                mapped[met_id] = mapping_details

            except:
                print ('erreur '+mol_chebi)

        mol_to_metabolites[mol_name] = mapped

    return mol_to_metabolites




def keep_closest_metabolites_only(mapping):
    """
    Kepps only mapped metabolites with the closest distance
    be careful when using this function with different kind of mappings merged together,
    it will keep the lowest value of distance, but it is hard to compare ditance between chebi mapping and lipidmaps mapping
    Args:
        arg1 (dict): mapping dictionnary
    Returns:
        dict: mapping dictionary with only mapped metabolites with smallest distance
    """
    clean_mapping ={}
    for mol_name in mapping:
        clean_mapping[mol_name] = {}
        distance = 99
        for met_id in mapping[mol_name]:
            new_distance = abs(mapping[mol_name][met_id]["distance"])
            if new_distance < distance:
                distance = new_distance

        for met_id in mapping[mol_name]:
            if abs(mapping[mol_name][met_id]["distance"]) <= distance:
                clean_mapping[mol_name][met_id] = mapping[mol_name][met_id]
            
    return clean_mapping


def merge_mappings(mapping_1, mapping_2):
    """
    Merge 2 mapping, adding their results together
    Args:
        arg1 (dict): mapping 1
        arg2 (dict): mapping 2
    Returns:
        dict: merge of both mappings with all found identifiers
    """
    if mapping_1 == {}:
        return mapping_2
    if mapping_2 == {}:
        return mapping_1
    for mol_name in mapping_2:
        for met_id in mapping_2[mol_name]:
            if met_id not in mapping_1[mol_name]:
                mapping_1[mol_name][met_id] = mapping_2[mol_name][met_id]
    return mapping_1