# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 14:33:03 2018

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

import argparse
import Parse as parse
import outputs as outputs
import Mapping as mapping
import csv
import collections
import json

#-------------------
#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("file_type", help="specify file type json for json files and tsv for tabulation separated files")
parser.add_argument("metabolomics_path", help= "Path to the metabolomics (including lipidomics) dataset")
parser.add_argument("network_metabolites_path", help = "Path to the network metabolites file")
parser.add_argument("output_path", help= "path to create the output file")
parser.add_argument("conf_file_path", help= "path to the conf file with the aliases")
parser.add_argument("mapping_types", help= "mappings that you wish to perfom. 1: exact multimapping, 2: chebi class mapping")
args = parser.parse_args()

#-------------------
#init arguments
#Type of files in which metabolomics data, network and output are provide
fileType = args.file_type
#aliases for labels
alias_dict = parse.parse_conf(args.conf_file_path)
#file containing the dataset
metabolomicsFile=args.metabolomics_path


#if the tsv option is selected, then we create a json file based on the tsv one
if fileType == 'tsv':
    # aliases
    OrderedDict = collections.OrderedDict
    src=metabolomicsFile
    dst=metabolomicsFile+'.json'
    header = list(alias_dict.values())
    data = []
    with open(src, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t', quotechar='"')
        next(reader)
        for row in reader:
            if row[0].strip()[0] == '#':  #
                continue
            row = filter(None, row)
            data.append(OrderedDict(zip(header, row)))
    with open(dst, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)
    metabolomicsFile = dst
    
    # ------------------------ START OF MY MODIFICATION OF THE CODE ---------------------------------------

    #get the metabolites in the network
    metabolites = parse.parse_metabolites_tsv(args.network_metabolites_path)
else:
    #get the metabolites in the network
    metabolites = parse.parse_metabolites_json(args.network_metabolites_path)
    # ------------------------ END OF MY MODIFICATION OF THE CODE -----------------------------------------


#dataset is produced based on the json file (provided by user or created based on txt file)
dataset = parse.parse_metabolomics_data_json(metabolomicsFile, alias_dict)

#------------------------- START ORIGINAL CODE ---------------------------------------------------------
""" 
#get the metabolites in the network
metabolites = parse.parse_metabolites_json(args.network_metabolites_path)
"""
# ------------------------- END ORIGINAL CODE ---------------------------------------------------------

### *** This returns only `None` for each metabolite it looked for in Chebi!!!!
metabolites_chebi_entities = mapping.convert_to_chebi_entities_once(metabolites)

#initialising empty dictionaries so that they exist when we try to merge the results in the end
multimapping = {}
chebi_mapping = {}

#determine mapping types and perform the different mappings
mapping_types = (args.mapping_types).split(",")
if "1" in mapping_types:
    multimapping = mapping.multimapping(dataset, metabolites)
if "2" in mapping_types:
    chebi_mapping = mapping.chebi_class_mapping(dataset, metabolites_chebi_entities)


final_merge = mapping.merge_mappings(multimapping, chebi_mapping)

clean_mapping = mapping.keep_closest_metabolites_only(final_merge)

#Create outPut files
if fileType == 'tsv':
    outputs.mapping_to_tsv(clean_mapping,args.output_path)
else:
    outputs.mapping_to_json(metabolomicsFile,metabolites, clean_mapping, args.output_path, alias_dict)