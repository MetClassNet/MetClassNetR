# -*- coding: utf-8 -*-
"""
Created on Oct 08 2019

@author: Florence Vinson and Fabien Jourdan

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
#coding: "utf8"
import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
exampledir=parentdir+'/example_data'
sys.path.insert(0,parentdir)

import Parse as parse
import Mapping as mapping

class Tests:
    # Set configuration file. To be tested
    alias_dict = parse.parse_conf(parentdir + '/conf.txt')
    # tests are performed on metabolites_4311_DB.json corresponding to Recon2.2
    metabolites = parse.parse_metabolites_json(exampledir + '/metabolites_4311_DB.json')
    metabolites_chebi_entities = mapping.convert_to_chebi_entities_once(metabolites)

    print('testing metabolite 15756')
    lipids = parse.parse_metabolomics_data_json(exampledir + '/data_15756.json', alias_dict)
    # initialising empty dictionaries so that they exist when we try to merge the results in the end
    multimapping = {}
    chebi_mapping = {}
    # Perform the two kinds of mapping
    multimapping = mapping.multimapping(lipids, metabolites)
    chebi_mapping = mapping.chebi_class_mapping(lipids, metabolites_chebi_entities)
    final_merge = mapping.merge_mappings(multimapping, chebi_mapping)
    clean_mapping = mapping.keep_closest_metabolites_only(final_merge)
    #check there is only one metabolite
    assert len(clean_mapping) == 1
    #check mapping retrieve 6 elements in the network
    assert len(clean_mapping['M1']) == 6
    #expected ones are: M_hdca_x;M_hdca_e;M_hdca_c;M_hdca_l;M_hdca_r;M_hdca_b
    assert not clean_mapping['M1']['M_hdca_x'] == None
    assert clean_mapping['M1']['M_hdca_x']['distance'] == -0.1
    assert clean_mapping['M1']['M_hdca_x']['mapping type'] == 'chebi class mapping'
    assert '7896' in clean_mapping['M1']['M_hdca_x']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_hdca_e'] == None
    assert clean_mapping['M1']['M_hdca_e']['distance'] == -0.1
    assert clean_mapping['M1']['M_hdca_e']['mapping type'] == 'chebi class mapping'
    assert '7896' in clean_mapping['M1']['M_hdca_e']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_hdca_c'] == None
    assert clean_mapping['M1']['M_hdca_c']['distance'] == -0.1
    assert clean_mapping['M1']['M_hdca_c']['mapping type'] == 'chebi class mapping'
    assert '7896' in clean_mapping['M1']['M_hdca_c']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_hdca_l'] == None
    assert clean_mapping['M1']['M_hdca_l']['distance'] == -0.1
    assert clean_mapping['M1']['M_hdca_l']['mapping type'] == 'chebi class mapping'
    assert '7896' in clean_mapping['M1']['M_hdca_l']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_hdca_r'] == None
    assert clean_mapping['M1']['M_hdca_r']['distance'] == -0.1
    assert clean_mapping['M1']['M_hdca_r']['mapping type'] == 'chebi class mapping'
    assert '7896' in clean_mapping['M1']['M_hdca_r']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_hdca_b'] == None
    assert clean_mapping['M1']['M_hdca_b']['distance'] == -0.1
    assert clean_mapping['M1']['M_hdca_b']['mapping type'] == 'chebi class mapping'
    assert '7896' in clean_mapping['M1']['M_hdca_b']['identifiers']['chebi']
    print('--test ok')

    print('testing metabolite 90488')
    lipids = parse.parse_metabolomics_data_json(exampledir + '/data_90488.json', alias_dict)
    # initialising empty dictionaries so that they exist when we try to merge the results in the end
    multimapping = {}
    chebi_mapping = {}
    # Perform the two kinds of mapping
    multimapping = mapping.multimapping(lipids, metabolites)
    chebi_mapping = mapping.chebi_class_mapping(lipids, metabolites_chebi_entities)
    final_merge = mapping.merge_mappings(multimapping, chebi_mapping)
    clean_mapping = mapping.keep_closest_metabolites_only(final_merge)
    #check there is only one metabolite
    assert len(clean_mapping) == 1
    #check mapping retrieve 6 elements in the network
    assert len(clean_mapping['M1']) == 4
    #expected ones are: M_pail_hs_c;M_pail_hs_g;M_pail_hs_n;M_pail_hs_r
    assert not clean_mapping['M1']['M_pail_hs_c'] == None
    assert clean_mapping['M1']['M_pail_hs_c']['distance'] == 1.1
    assert clean_mapping['M1']['M_pail_hs_c']['path'] == ['90488', '16749']
    assert clean_mapping['M1']['M_pail_hs_c']['mapping type'] == 'chebi class mapping'
    assert '57880' in clean_mapping['M1']['M_pail_hs_c']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_pail_hs_g'] == None
    assert clean_mapping['M1']['M_pail_hs_g']['distance'] == 1.1
    assert clean_mapping['M1']['M_pail_hs_g']['path'] == ['90488', '16749']
    assert clean_mapping['M1']['M_pail_hs_g']['mapping type'] == 'chebi class mapping'
    assert '57880' in clean_mapping['M1']['M_pail_hs_g']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_pail_hs_n'] == None
    assert clean_mapping['M1']['M_pail_hs_n']['distance'] == 1.1
    assert clean_mapping['M1']['M_pail_hs_n']['path'] == ['90488', '16749']
    assert clean_mapping['M1']['M_pail_hs_n']['mapping type'] == 'chebi class mapping'
    assert '57880' in clean_mapping['M1']['M_pail_hs_n']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_pail_hs_r'] == None
    assert clean_mapping['M1']['M_pail_hs_r']['distance'] == 1.1
    assert clean_mapping['M1']['M_pail_hs_r']['path'] == ['90488', '16749']
    assert clean_mapping['M1']['M_pail_hs_r']['mapping type'] == 'chebi class mapping'
    assert '57880' in clean_mapping['M1']['M_pail_hs_r']['identifiers']['chebi']
    print('--test ok')

    print('testing metabolite 36023')
    lipids = parse.parse_metabolomics_data_json(exampledir + '/data_36023.json', alias_dict)
    # initialising empty dictionaries so that they exist when we try to merge the results in the end
    multimapping = {}
    chebi_mapping = {}
    # Perform the two kinds of mapping
    multimapping = mapping.multimapping(lipids, metabolites)
    chebi_mapping = mapping.chebi_class_mapping(lipids, metabolites_chebi_entities)
    final_merge = mapping.merge_mappings(multimapping, chebi_mapping)
    clean_mapping = mapping.keep_closest_metabolites_only(final_merge)
    #check there is only one metabolite
    assert len(clean_mapping) == 1
    #check mapping retrieve 6 elements in the network
    assert len(clean_mapping['M1']) == 3
    #expected ones are: M_vacc_e;M_vacc_b;M_vacc_c
    assert not clean_mapping['M1']['M_vacc_e'] == None
    assert clean_mapping['M1']['M_vacc_e']['distance'] == -1.1
    assert clean_mapping['M1']['M_vacc_e']['path'] == ['36023', '28727']
    assert clean_mapping['M1']['M_vacc_e']['mapping type'] == 'chebi class mapping'
    assert '30828' in clean_mapping['M1']['M_vacc_e']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_vacc_e'] == None
    assert clean_mapping['M1']['M_vacc_e']['distance'] == -1.1
    assert clean_mapping['M1']['M_vacc_e']['path'] == ['36023', '28727']
    assert clean_mapping['M1']['M_vacc_e']['mapping type'] == 'chebi class mapping'
    assert '30828' in clean_mapping['M1']['M_vacc_e']['identifiers']['chebi']
    assert not clean_mapping['M1']['M_vacc_e'] == None
    assert clean_mapping['M1']['M_vacc_e']['distance'] == -1.1
    assert clean_mapping['M1']['M_vacc_e']['path'] == ['36023', '28727']
    assert clean_mapping['M1']['M_vacc_e']['mapping type'] == 'chebi class mapping'
    assert '30828' in clean_mapping['M1']['M_vacc_e']['identifiers']['chebi']

    print('--test ok')

