# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 09:21:30 2018

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

# ---------------------- START OF MY MODIFICATIONS OF THE CODE --------------------
import pandas as pd
import re
"""
===============================================================================
==================== Functions to parse the TSV files =========================
===============================================================================
"""
def parse_metabolites_tsv(metabolite_table_path):
    """
    Parses a TSV file containing all the metabolites and their attributes. The file
    must contain 4 columns: id, name, formula, and chebi
    
    Args:
        arg1 (str): Path to the metabolite TSV file
    
    Returns: 
        dict: metabolite dictionary {networkid : dictionary of identifiers}
    """
    metabolites = {}
    metaboliteTable = pd.read_csv(metabolite_table_path, sep = "\t")
    metaboliteTable.columns = metaboliteTable.columns.str.lower()
    chebis = metaboliteTable.loc[:, "chebi"]

    # chebi idetifiers must have only the number, not the "CHEBI: " part 
    for row in list(range(0, len(chebis))):
        # get chebi
        c = chebis[row]

        # split string (it is needed when there are multpile chebi IDs)
        ids = c.split(";")

        cleanChebi = []

        for i in list(range(0, len(ids))):
            id = ids[i]

            # get chebi number
            number = re.search(r".*?(\d+)", id).group(1)

            # add chebi number
            cleanChebi.append(str(number))

        metabolites[metaboliteTable.loc[row, "id"]] = { "name": [metaboliteTable.loc[row, "name"]], "formula": [metaboliteTable.loc[row, "formula"]], "chebi" : cleanChebi } 

    return metabolites 
# ---------------------- END OF MY MODIFICATIONS OF THE CODE --------------------


import json
#from rdkit import Chem


"""
===============================================================================
==================== Functions to parse the json files ========================
===============================================================================
"""

def parse_metabolites_json(metabolite_json_path):
    """
    Parses a json file containing all the metabolites and their different identifiers as such
    
    {
    "success":true,
    "results":[
                {
                "id":"106010",
                 "name":"10,11-dihydro-12R-hydroxy-leukotriene C4",
                 "dbIdentifier":"M_CE4993_x","chemicalFormula":"C30H47N3O10S",
                 "idin":[
                         {
                         "extDBName":"SBO",
                         "extID":"SBO:0000247",
                         "origin":"SBML File",
                         "score":"1"
                         }
                         
                        ]
                 }
                 
             ]
    }
    
    Check that all the text part in the files are between quote, or the json library won't be able to parse the file
    
    Args:
        arg1 (str): Path to the metabolite json file.
    
    Returns: 
        dict: metabolite dictionary {networkid : dictionary of identifiers}
    """
    metabolites = {}
    with open(metabolite_json_path, "r") as metabolites_file:
        data = json.load(metabolites_file)
        met_list = data["results"]
        for metabolite in met_list:
            met_networkid = metabolite["dbIdentifier"]
            
            met_name = [metabolite["name"]]
            met_formula = [metabolite["chemicalFormula"]]
            met_chebi = []
            met_kegg = []
            met_pubchem = []
            met_inchi = []
            met_hmdb = []
            met_smiles = []
            met_lmids = []
            met_swissids = []
            iddb = [metabolite["id"]]
            
            
            for identifier in metabolite["idin"]:
                
                if identifier["extDBName"] == "chebi":
                    chebi = identifier["extID"].replace(" ", "").replace("\n","")
                    if "chebi" in chebi.lower():
                        chebi = chebi[6:]
                    met_chebi.append(chebi)

                elif identifier["extDBName"] == "kegg.compound":
                    met_kegg.append(identifier["extID"].replace(" ", "").replace("\n",""))
                    
                elif identifier["extDBName"] == "pubchem.compound":
                    met_pubchem.append(identifier["extID"].replace(" ", "").replace("\n",""))
                
                elif identifier["extDBName"] == "inchi":
                    inchi = identifier["extID"].replace(" ", "").replace("\n","")
                    if "inchi=" in inchi.lower():
                        inchi = inchi[6:]
#                    if inchi[-4:] == "/p-1":
#                       inchi = inchi[:-4]
                    if "1S/" in inchi:
                        inchi = inchi[3:]
                    elif "1/" in met_inchi:
                        inchi = inchi[2:]
                    
                    met_inchi.append(inchi)
                        

                elif identifier["extDBName"] == "hmdb":
                    met_hmdb.append(identifier["extID"].replace(" ", "").replace("\n",""))                 
                
                elif identifier["extDBName"] == "smiles":
                    smiles = identifier["extID"].replace(" ", "").replace("\n","")
                    if "smiles:" in smiles.lower():
                        smiles = smiles[7:]
#                    try:  
#                        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
#                    except: smiles = smiles
                    met_smiles.append(smiles)
                    
                
                elif identifier["extDBName"] == "lipidmaps":
                    met_lmids.append(identifier["extID"].replace(" ", "").replace("\n","")) 
                
                elif identifier["extDBName"] == "swisslipids":
                    met_swissids.append(identifier["extID"].replace(" ", "").replace("\n",""))
            
            if "idout" in metabolites:
                for identifier in metabolites["idout"]:
                
                    if identifier["extDBName"] == "chebi":
                        chebi = identifier["extID"].replace(" ", "").replace("\n","")
                        if "chebi" in chebi.lower():
                            chebi = chebi[6:]
                        met_chebi.append(chebi)
    
                    elif identifier["extDBName"] == "kegg.compound":
                        met_kegg.append(identifier["extID"].replace(" ", "").replace("\n",""))
                        
                    elif identifier["extDBName"] == "pubchem.compound":
                        met_pubchem.append(identifier["extID"].replace(" ", "").replace("\n",""))
                    
                    elif identifier["extDBName"] == "inchi":
                        inchi = identifier["extID"].replace(" ", "").replace("\n","")
                        if "inchi=" in inchi.lower():
                            inchi = inchi[6:]
                        if inchi[-4:] == "/p-1":
                            inchi = inchi[:-4]
                        if "1S/" in inchi:
                            inchi = inchi[3:]
                        elif "1/" in met_inchi:
                            inchi = inchi[2:]
                        
                        met_inchi.append(inchi)
                            
    
                    elif identifier["extDBName"] == "hmdb":
                        met_hmdb.append(identifier["extID"].replace(" ", "").replace("\n",""))                 
                    
                    elif identifier["extDBName"] == "smiles":
                        smiles = identifier["extID"].replace(" ", "").replace("\n","")
                        if "smiles:" in smiles.lower():
                            smiles = smiles[7:]
    #                    try:  
    #                        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    #                    except: smiles = smiles
                        met_smiles.append(smiles)
                        
                    
                    elif identifier["extDBName"] == "lipidmaps":
                        met_lmids.append(identifier["extID"].replace(" ", "").replace("\n","")) 
                    
                    elif identifier["extDBName"] == "swisslipids":
                        met_swissids.append(identifier["extID"].replace(" ", "").replace("\n","")) 
                
            metabolites[met_networkid] = {"name" : met_name, "formula" : met_formula, "chebi" : met_chebi, 
                                        "kegg" : met_kegg, "pubchem" :met_pubchem, "inchi" : met_inchi, 
                                        "hmdb" : met_hmdb, "smiles": met_smiles, "lmid": met_lmids, "swissid": met_swissids, "iddb" : iddb}
    
    return metabolites 





def parse_conf(conf_file_path):
    """
        Parsing the conf file to determine the keys that are used in the json file.
    Args:
        arg1 (str): path to the configuration file for json
    Returns:
        dict: dictionary of key aliases that will be used to parse the network metabolites and dataset files
    """
    alias_dict = {}
    with open(conf_file_path, "r") as conf:
        for line in conf:
            if "name alias:" in line:
                name_alias = line[12:]
                if name_alias[0] == " ":
                    name_alias = name_alias[1:]
                alias_dict["name_alias"]=name_alias.replace("\n","")

                
            elif "formula alias:" in line:
                formula_alias = line[14:]
                if formula_alias[0]== " ":
                    formula_alias = formula_alias[1:]
                alias_dict["formula_alias"]=formula_alias.replace("\n","")

                    
            elif "chebi alias:" in line:
                chebi_alias = line[12:]
                if chebi_alias[0] == " ":
                    chebi_alias = chebi_alias[1:]
                alias_dict["chebi_alias"]=chebi_alias.replace("\n","")
            
            elif "pubchem alias:" in line:
                pubchem_alias = line[14:]
                if pubchem_alias[0] == " ":
                    pubchem_alias = pubchem_alias[1:]
                alias_dict["pubchem_alias"]=pubchem_alias.replace("\n","")
                    
            elif "kegg alias:" in line:
                kegg_alias = line[11:]
                if kegg_alias[0] == " ":
                    kegg_alias = kegg_alias[1:]
                alias_dict["kegg_alias"]=kegg_alias.replace("\n","")

    
            elif "inchi alias:" in line:
                inchi_alias = line[12:]
                if inchi_alias[0] == " ":
                    inchi_alias = inchi_alias[1:]
                alias_dict["inchi_alias"]=inchi_alias.replace("\n","")

    
            elif "hmdb alias:" in line:
                hmdb_alias = line[11:]
                if hmdb_alias[0] == " ":
                    hmdb_alias = hmdb_alias[1:]
                alias_dict["hmdb_alias"]=hmdb_alias.replace("\n","")

            
            elif "smiles alias:" in line:
                smiles_alias = line[13:]
                if smiles_alias[0] == " ":
                    smiles_alias = smiles_alias[1:]
                alias_dict["smiles_alias"]=smiles_alias.replace("\n","")

            
            elif "lipidmaps alias:" in line:
                lmid_alias = line[16:]
                if lmid_alias[0] == " ":
                    lmid_alias = lmid_alias[1:]
                alias_dict["lipidmaps_alias"]=lmid_alias.replace("\n","")

    
            elif "swisslipids alias:" in line:
                swissid_alias = line[18:]
                if swissid_alias[0] == " ":
                    swissid_alias = swissid_alias[1:]
                alias_dict["swisslipids_alias"]=swissid_alias.replace("\n","")

    return alias_dict



def parse_metabolomics_data_json(metabolomics_data_json_path, alias_dict):
    """
    Parses a json file of metabolites obtained in metabolomics, json key are determined by using the conf_dict obtained by parsing the conf_file
       
    Args:
        arg1 (str): path to json file of metabolites
        arg2 (dict): configuration dictionary obtained using parse_conf_json
    Returns:
        dict: Dictionary: {metabolite_name : dictionary of identifiers}

    """
    name_alias = alias_dict["name_alias"]
    formula_alias = alias_dict["formula_alias"]
    chebi_alias = alias_dict["chebi_alias"]
    kegg_alias = alias_dict["kegg_alias"]
    pubchem_alias = alias_dict["pubchem_alias"]
    hmdb_alias = alias_dict["hmdb_alias"]
    inchi_alias = alias_dict["inchi_alias"]
    smiles_alias = alias_dict["smiles_alias"]
    lmid_alias = alias_dict["lipidmaps_alias"]
    swissid_alias = alias_dict["swisslipids_alias"]
    
    
    metabolites = {}
    with open(metabolomics_data_json_path, "r") as metabolomics_file:
        data = json.load(metabolomics_file)

        for metabolite in data:
            
            metabolite_name = metabolite[name_alias]
            
            try: metabolite_formula = metabolite[formula_alias].replace(" ","").replace("\n","")
            except: metabolite_formula = ""
            
            try:
                metabolite_chebi = metabolite[chebi_alias].replace(" ","").replace("\n","")
                if "chebi" in metabolite_chebi.lower():
                    metabolite_chebi = metabolite_chebi[6:]
            except: metabolite_chebi = ""
            
            try:metabolite_kegg = metabolite[kegg_alias].replace(" ","").replace("\n","")
            except: metabolite_kegg = ""
            
            try: metabolite_pubchem = metabolite[pubchem_alias].replace(" ","").replace("\n","")
            except: metabolite_pubchem = ""
            
            try: metabolite_hmdb = metabolite[hmdb_alias].replace(" ","").replace("\n","")
            except: metabolite_hmdb = ""
            
            try :
                metabolite_inchi = metabolite[inchi_alias].replace(" ","").replace("\n","")
                if "inchi=" in metabolite_inchi.lower():
                    metabolite_inchi = metabolite_inchi[6:]
                if metabolite_inchi[-4:] == "/p-1":
                    metabolite_inchi = metabolite_inchi[:-4]
                if "1S/" in metabolite_inchi:
                    metabolite_inchi = metabolite_inchi[3:]
                elif "1/" in metabolite_inchi:
                    metabolite_inchi = metabolite_inchi[2:]
                
            except: metabolite_inchi = ""
            
            try :
                metabolite_smiles = metabolite[smiles_alias].replace(" ","").replace("\n","")
                if "smiles:" in metabolite_inchi.lower():
                    metabolite_smiles = metabolite_smiles[7:]
#                    try:
#                        metabolite_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(metabolite_smiles))
#                    except: metabolite_smiles = metabolite_smiles
            except: metabolite_smiles = ""
                    
            try: metabolite_lmid = metabolite[lmid_alias].replace(" ","").replace("\n","")
            except: metabolite_lmid = ""
            
            try: metabolite_swissid = metabolite[swissid_alias].replace(" ","").replace("\n","")
            except: metabolite_swissid = ""
            




            metabolites[metabolite_name] ={"name" : metabolite_name, "formula" : metabolite_formula, "chebi" : metabolite_chebi,
                                    "kegg" : metabolite_kegg, "pubchem" :metabolite_pubchem, "inchi" : metabolite_inchi,
                                    "hmdb" : metabolite_hmdb, "smiles": metabolite_smiles, "lmid" : metabolite_lmid, "swissid": metabolite_swissid}
            
        
    return metabolites


"""
===============================================================================
==================== Functions to parse the tsv files =========================
===============================================================================
"""


def parse_conf_tsv(conf_file_path):
    """
    Parsing the conf file to determine the names of the columns.
    Args:
        arg1 (str): path to the configuration file for tsv
    Returns:
        dict: dictionary of column names that will be used to parse the network metabolites and dataset files
    """
    met_conf_dict = {}
    dataset_conf_dict = {}
    with open(conf_file_path, "r") as conf:
        for line in conf:
            
            if "met multiple identifiers separator" in line:
                sep = line[35:]
                if sep[0] == " ":
                    sep = sep[1:]
                met_conf_dict["sep"] = sep.replace("\n","")
            
            #First read the part to determine the columns in the metabolite file
            elif "met id column" in line:
               met_id_col = line[14:]
               if met_id_col[0] == " ":
                   met_id_col = met_id_col[1:]
               met_conf_dict["met_id_col"] = met_id_col.replace("\n","")         
            
            elif "met name column" in line:
               met_name_col = line[16:]
               if met_name_col[0] == " ":
                   met_name_col = met_name_col[1:]
               met_conf_dict["met_name_col"] = met_name_col.replace("\n","")
               
            elif "met formula column" in line:
               met_formula_col = line[19:]
               if met_formula_col[0] == " ":
                   met_formula_col = met_formula_col[1:]
               met_conf_dict["met_formula_col"] = met_formula_col.replace("\n","") 
               
            elif "met chebi column" in line :
               met_chebi_col = line[17:]
               if met_chebi_col[0] == " ":
                   met_chebi_col = met_chebi_col[1:]
               met_conf_dict["met_chebi_col"] = met_chebi_col.replace("\n","") 
                   
            elif "met kegg column" in line:
               met_kegg_col = line[16:]
               if met_kegg_col[0] == " ":
                   met_kegg_col= met_kegg_col[1:]
               met_conf_dict["met_kegg_col"] = met_kegg_col.replace("\n","") 
                   
            elif "met pubchem column" in line:
               met_pubchem_col = line[19:]
               if met_pubchem_col[0] == " ":
                   met_pubchem_col = met_pubchem_col[1:]
               met_conf_dict["met_pubchem_col"] = met_pubchem_col.replace("\n","") 
                  
            elif "met inchi column" in line:
               met_inchi_col = line[17:]
               if met_inchi_col[0] == " ":
                   met_inchi_col = met_inchi_col[1:]
               met_conf_dict["met_inchi_col"] = met_inchi_col.replace("\n","") 
                   
            elif "met hmdb column" in line:
               met_hmdb_col = line[16:]
               if met_hmdb_col[0]== " ":
                   met_hmdb_col = met_hmdb_col[1:]
               met_conf_dict["met_hmdb_col"] = met_hmdb_col.replace("\n","")
                   
            elif "met smiles column" in line:
               met_smiles_col = line[18:]
               if met_smiles_col[0] == " ":
                   met_smiles_col = met_smiles_col [1:]
               met_conf_dict["met_smiles_col"] = met_smiles_col.replace("\n","") 

            elif "met lipidmaps column" in line:
               met_lmid_col = line[21:]
               if met_lmid_col[0] == " ":
                   met_lmid_col = met_lmid_col [1:]
               met_conf_dict["met_lmid_col"] = met_lmid_col.replace("\n","")

            elif "met swisslipids column" in line:
               met_swissid_col = line[23:]
               if met_swissid_col[0] == " ":
                   met_swissid_col = met_swissid_col[1:]
               met_conf_dict["met_swissid_col"] = met_swissid_col.replace("\n","")

                 
            #Now we determine the columns names for the dataset file
            
            elif "dataset name column" in line:
               dataset_name_col = line[18:]
               if dataset_name_col[0] == " ":
                   dataset_name_col = dataset_name_col[1:]
               #print(dataset_name_col.replace("\n",""))
               dataset_conf_dict["dataset_name_col"] = dataset_name_col.replace("\n","")
               
            elif "dataset formula column" in line:
               dataset_formula_col = line[21:]
               if dataset_formula_col[0] == " ":
                   dataset_formula_col = dataset_formula_col[1:]
               dataset_conf_dict["dataset_formula_col"] = dataset_formula_col.replace("\n","")
               
            elif "dataset chebi column" in line :
               dataset_chebi_col = line[19:]
               if dataset_chebi_col[0] == " ":
                   dataset_chebi_col = dataset_chebi_col[1:]
               dataset_conf_dict["dataset_chebi_col"] = dataset_chebi_col.replace("\n","")
                   
            elif "dataset kegg column" in line:
               dataset_kegg_col = line[18:]
               if dataset_kegg_col[0] == " ":
                   dataset_kegg_col= dataset_kegg_col[1:]
               dataset_conf_dict["dataset_kegg_col"] = dataset_kegg_col.replace("\n","")
                   
            elif "dataset pubchem column" in line:
               dataset_pubchem_col = line[21:]
               if dataset_pubchem_col[0] == " ":
                   dataset_pubchem_col = dataset_pubchem_col[1:]
               dataset_conf_dict["dataset_pubchem_col"] = dataset_pubchem_col.replace("\n","")
                  
            elif "dataset inchi column" in line:
               dataset_inchi_col = line[19:]
               if dataset_inchi_col[0] == " ":
                   dataset_inchi_col = dataset_inchi_col[1:]
               dataset_conf_dict["dataset_inchi_col"] = dataset_inchi_col.replace("\n","")
                   
            elif "dataset hmdb column" in line:
               dataset_hmdb_col = line[18:]
               if dataset_hmdb_col[0]== " ":
                   dataset_hmdb_col = dataset_hmdb_col[1:]
               dataset_conf_dict["dataset_hmdb_col"] = dataset_hmdb_col.replace("\n","")
                   
            elif "dataset smiles column" in line:
               dataset_smiles_col = line[20:]
               if dataset_smiles_col[0] == " ":
                   dataset_smiles_col = dataset_smiles_col [1:]
               dataset_conf_dict["dataset_smiles_col"] = dataset_smiles_col.replace("\n","")

            elif "dataset lipidmaps column" in line:
               dataset_lmid_col = line[23:]
               if dataset_lmid_col[0] == " ":
                   dataset_lmid_col = dataset_lmid_col [1:]
               dataset_conf_dict["dataset_lmid_col"] = dataset_lmid_col.replace("\n","")

            elif "dataset swisslipids column" in line:
               dataset_swissid_col = line[25:]
               if dataset_swissid_col[0] == " ":
                   dataset_swissid_col = dataset_swissid_col[1:]
               dataset_conf_dict["dataset_swissid_col"] = dataset_swissid_col.replace("\n","")
    
    
    return met_conf_dict, dataset_conf_dict




#-------Deprecated
def parse_networkMetabolites_tsv(metabolites_tsv_path, conf_dict):
    
    """
    Args:
        arg1 (str): Path to the metabolite tabulated file.
        arg2 (dict): Dictionnary of the names of the columns obtained from the conf file    
    Returns: 
        dict: Metabolite dictionary {networkid : dictionary of identifiers}
    """

    metabolites = {}
    with open(metabolites_tsv_path, mode="r") as metabolites_file:
        header_line = next(metabolites_file)
        headers = header_line[:-1].split("\t")
        #setting the column indexes, can be done using a conf file
        sep = conf_dict["sep"]
        #poor hacking by Fabien to allow the code to work
        global id_col
        id_col = 1
        for i in range(len(headers)):
            if headers[i] == conf_dict["met_id_col"]:
                id_col = i
            if headers[i] == conf_dict["met_name_col"] : 
                name_col = i
            if headers[i] == conf_dict["met_formula_col"]:
                formula_col = i
            if headers[i] == conf_dict["met_chebi_col"]:
                chebi_col = i
            if headers[i] == conf_dict["met_kegg_col"]:
                kegg_col = i
            if headers[i] == conf_dict["met_pubchem_col"]:
                pubchem_col = i
            if headers[i] == conf_dict["met_inchi_col"]:
                inchi_col = i
            if headers[i] == conf_dict["met_hmdb_col"]:
                hmdb_col = i
            if headers[i] == conf_dict["met_smiles_col"]:
                smiles_col = i
            if headers[i] == conf_dict["met_lmid_col"]:
                lmid_col = i
            if headers[i] == conf_dict["met_swissid_col"]:
                swissid_col = i

            
        for line in metabolites_file:
            met_data = line[:-1].split("\t") #removing \n last character of each line and splitting with tabulations
            met_networkid = met_data[id_col]
            met_names = []
            met_formulas = []
            met_chebis = []
            met_keggs = []
            met_pubchems = []
            met_inchis = []
            met_hmdbs = []
            met_smiles = []
            met_lmids = []
            met_swissids = []
            
            
            try:
                names = met_data[name_col]
                for name in names.split(sep):
                    name = name.replace("\n","")
                    if name != "" :
                        met_names.append(name)
            except: pass
                    
            try:
                formulas = met_data[formula_col]
                for formula in formulas.split(sep):
                    formula = formula.replace(" ", "").replace("\n","")
                    if formula != "":
                        met_formulas.append(formula)
            except: pass
            
            try:
                chebis = met_data[chebi_col]
                for chebi in chebis.split(sep):
                    chebi = chebi.replace(" ", "").replace("\n","")
                    if "chebi:" in chebi.lower():
                        chebi = chebi[6:]
                    if chebi != "":
                        while chebi[-1] not in "0123456789":
                            chebi = chebi[:-1]
                        met_chebis.append(chebi)
            except: pass
          
                    
            try:
                keggs = met_data[kegg_col]
                for kegg in keggs.split(sep):
                    kegg = kegg.replace(" ", "").replace("\n","")
                    if kegg != "":
                        met_keggs.append(kegg)
            except: pass
                      
            try:
                pubchems = met_data[pubchem_col]
                for pubchem in pubchems.split(sep):
                    pubchem= pubchem.replace(" ", "").replace("\n","")
                    if pubchem != "":
                        met_pubchems.append(pubchem)
            except: pass
                          
            try:
                inchis = met_data[inchi_col]
                for inchi in inchis.split(sep):
                    inchi = inchi.replace(" ", "").replace("\n","")
                    if "inchi=" in inchi.lower():
                        inchi = inchi[6:]
                    if inchi[-4:] == "/p-1": #/p-1 is the conjugate base, I remove them so I don't miss a match just because its the base instead of the acid (or reverse)
                        inchi = inchi[:-4]
                    if "1S/" in inchi:
                        inchi = inchi[3:]
                    elif "1/" in inchi:
                        inchi = inchi[2:]
                    if inchi != "":
                        met_inchis.append(inchi)
            except: pass
                                  
                                    
            try:
                hmdbs = met_data[hmdb_col].split(sep)
                for hmdb in hmdbs:
                    hmdb = hmdb.replace(" ", "").replace("\n","")
                    if hmdb != "":
                        met_hmdbs.append(hmdb)
            except: pass
            
            try:
                smiles = met_data[smiles_col].split(sep)
                for smile in smiles:
                    smile = smile.replace(" ", "").replace("\n","")
                if "smiles:" in smile.lower():
                    smile = smile[7:]
#               try:
#                   smile = Chem.MolToSmiles(Chem.MolFromSmiles(smile)) 
#               except: smile = smile
                if smile != "":
                    met_smiles.append(smile)
            except: pass
            
            try:
                lmids = met_data[lmid_col].split(sep)
                for lmid in lmids:
                    lmid = lmid.replace(" ", "").replace("\n","")
                    if lmid != "":
                        met_hmdbs.append(lmid)
            except: pass
                    
            try:
                swissids = met_data[swissid_col].split(sep)
                for swissid in swissids:
                    swissid = swissid.replace(" ", "").replace("\n","")
                    if swissid != "":
                        met_swissids.append(swissid)
            except: pass
            
            metabolites[met_networkid] = {"name" : met_names, "formula" : met_formulas, "chebi" : met_chebis, 
                                           "kegg" : met_keggs, "pubchem" :met_pubchems, "inchi" : met_inchis, 
                                           "hmdb" : met_hmdbs, "smiles": met_smiles, "lmid": met_lmids, "swissid": met_swissids}
            
            
            
    
    return metabolites



def parse_metabolomics_data_tsv(metabolomics_data_tsv_path, conf_dict):
    """
    Parses a tabulated file of metabolites (metabolomics data), columns are determined by using the conf_dict obtained by parsing the conf_file
       
    Args:
        arg1 (str): Tabulated file of dataset path
        arg2 (dict): Configuration dictionary obtained using parse_conf_tsv
    Returns:
        dict: Dictionary: {metabolite_name : dictionary of identifiers}

    """
    metabolites = {}
    with open(metabolomics_data_tsv_path, "r") as dataset_file:
        header_line = next(dataset_file)
        headers = header_line[:-1].split("\t")
        #setting the column indexes, can be done using a conf file        
        for i in range(len(headers)):
            if headers[i] == conf_dict["dataset_name_col"]:
                name_col = i
            if headers[i] == conf_dict["dataset_formula_col"]:
                formula_col = i
            if headers[i] == conf_dict["dataset_chebi_col"]:
                chebi_col = i
            if headers[i] == conf_dict["dataset_kegg_col"]:
                kegg_col = i
            if headers[i] == conf_dict["dataset_pubchem_col"]:
                pubchem_col = i
            if headers[i] == conf_dict["dataset_inchi_col"]:
                inchi_col = i
            if headers[i] == conf_dict["dataset_hmdb_col"]:
                hmdb_col = i
            if headers[i] == conf_dict["dataset_smiles_col"]:
                smiles_col = i
            if headers[i] == conf_dict["dataset_lmid_col"]:
                lmid_col = i
            if headers[i] == conf_dict["dataset_swissid_col"]:
                swissid_col = i
        for line in dataset_file:
    
            data = (line[:-1].split("\t"))

            try:
                metabolite_name = data[name_col]
            except:
                metabolite_name = ""
            try:
                metabolite_formula = data[formula_col]
            except: metabolite_formula = ""
            
            try:
                metabolite_chebi = data[chebi_col].replace(" ", "").replace("\n","")
                if metabolite_chebi != "":
                    if "chebi" in metabolite_chebi.lower():
                        metabolite_chebi = metabolite_chebi[6:]
                    while metabolite_chebi[-1] not in "0123456789": #used to remove weird character that happen to be at the end of some chebi strings in the data...
                        metabolite_chebi = metabolite_chebi[:-1]
            except: metabolite_chebi = ""
        
            try:
                metabolite_kegg = data[kegg_col].replace(" ", "").replace("\n","")
            except: metabolite_kegg = ""
        
            try:
                metabolite_pubchem = data[pubchem_col].replace(" ", "").replace("\n","")
            except: metabolite_pubchem = ""

            try:
                metabolite_inchi = data[inchi_col]
                metabolite_inchi = metabolite_inchi.replace(" ", "").replace("\n","")
                if "inchi=" in metabolite_inchi.lower():
                    metabolite_inchi = metabolite_inchi[6:]
                if metabolite_inchi[-4:] == "/p-1":
                    metabolite_inchi = metabolite_inchi[:-4]
                if "1S/" in metabolite_inchi:
                    metabolite_inchi = metabolite_inchi[3:]
                elif "1/" in metabolite_inchi:
                    metabolite_inchi = metabolite_inchi[2:]
            except: metabolite_inchi = ""
            
            try:
                metabolite_hmdb = data[hmdb_col].replace(" ", "").replace("\n","")
            except: metabolite_hmdb = ""
            
            try:
                metabolite_smiles = data[smiles_col].replace(" ", "").replace("\n","")
                if "smiles" in metabolite_smiles.lower():
                    metabolite_smiles = metabolite_smiles[6:]
    #                try:
    #                    metabolite_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(metabolite_smiles))
    #                except: metabolite_smiles = metabolite_smiles
            except: metabolite_smiles = ""

            try:
                metabolite_lmid = data[lmid_col].replace(" ", "").replace("\n","")
            except: metabolite_lmid = ""
            
            try:
                metabolite_swissid = data[swissid_col].replace(" ", "").replace("\n","")
                if "slm:" in metabolite_swissid.lower():
                    metabolite_swissid = metabolite_swissid[4:]
            except: metabolite_swissid = ""

            
            metabolites[metabolite_name] ={"name" : metabolite_name, "formula" : metabolite_formula, "chebi" : metabolite_chebi,
                                    "kegg" : metabolite_kegg, "pubchem" :metabolite_pubchem, "inchi" : metabolite_inchi,
                                    "hmdb" : metabolite_hmdb, "smiles": metabolite_smiles, "lmid" : metabolite_lmid, "swissid": metabolite_swissid}
            

    return metabolites

    