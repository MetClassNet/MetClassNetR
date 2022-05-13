# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 11:46:53 2018

@author: atmoreau

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

import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import Mapping as mapping
import unittest


class TestMapping(unittest.TestCase):
    """
    For each parsing function, tests are done comparing that the data parsed for a metabolite from the test files
    is conform to what is expected. 
    """
    
    lipids = {
                "lipid1" : {"name":"lipid1", "formula":"formula1","chebi":"00001","inchi":"1S/inchi1","hmdb":"hmdb1", 
                            "kegg":"kegg1", "pubchem":"pubchem1","lmid":"LMFA01010001", "swissid" : "swissid1", "smiles":"smiles1"},
                "lipid2" : {"name":"lipid2", "formula":"formula2","chebi":"00002","inchi":"1S/inchi2","hmdb":"hmdb2", 
                            "kegg":"kegg2", "pubchem":"pubchem2","lmid":"LMFA01010002", "swissid" : "swissid2", "smiles":"smiles2"},
                "lipid3" : {"name":"lipid3", "formula":"formula3","chebi":"16196","inchi":"1S/inchi3","hmdb":"hmdb3", 
                            "kegg":"kegg3", "pubchem":"pubchem3","lmid":"LMFA01010003", "swissid" : "swissid3", "smiles":"smiles3"},
                "lipid4" : {"name":"lipid4", "formula":"formula4","chebi":"00004","inchi":"1S/inchi4","hmdb":"hmdb4", 
                            "kegg":"kegg4", "pubchem":"pubchem4","lmid":"LMFA01010004", "swissid" : "swissid4", "smiles":"smiles4"},
                "lipid5" : {"name":"lipid5", "formula":"formula5","chebi":"00005","inchi":"1S/inchi5","hmdb":"hmdb5", 
                            "kegg":"kegg5", "pubchem":"pubchem5","lmid":"LMFA01010005", "swissid" : "swissid5", "smiles":"smiles5"},
                "lipid6" : {"name":"lipid6", "formula":"formula6","chebi":"00006","inchi":"1S/inchi6","hmdb":"hmdb6", 
                            "kegg":"kegg6", "pubchem":"pubchem6","lmid":"LMFA01010006", "swissid" : "swissid6", "smiles":"smiles6"}    
            }
    
    
    metabolites = {
                    "met1" : {"name":["lipid1"], "formula":["formula1"],"chebi":["00001"],"inchi":["1S/inchi1"],"hmdb":["hmdb1"], 
                                "kegg":["kegg1"], "pubchem":["pubchem1"],"lmid":["LMFA01010001"], "swissid" : ["swissid1"], "smiles":["smiles1"]},
                    "met2" : {"name":["lipid2"], "formula":["formula2"],"chebi":["00002"],"inchi":["1S/inchi2"],"hmdb":["hmdb2"], 
                                "kegg":["kegg2"], "pubchem":["pubchem2"],"lmid":["LMFA01010002"], "swissid" : ["swissid2"], "smiles":["smiles2"]},
                    "met3" : {"name":["lipid3"], "formula":["formula3"],"chebi":["35366"],"inchi":["1S/inchi3"],"hmdb":["hmdb3"], 
                                "kegg":["kegg3"], "pubchem":["pubchem3"],"lmid":["LMFA01010003"], "swissid" : ["swissid3"], "smiles":["smiles3"]},
                    "met4" : {"name":["lipid4"], "formula":["formula4"],"chebi":["00004"],"inchi":["1S/inchi4"],"hmdb":["hmdb4"], 
                                "kegg":["kegg4"], "pubchem":["pubchem4"],"lmid":["LMFA01010004"], "swissid" : ["swissid4"], "smiles":["smiles4"]},
                    "met5" : {"name":["lipid5"], "formula":["formula5"],"chebi":["00005"],"inchi":["1S/inchi5"],"hmdb":["hmdb5"], 
                                "kegg":["kegg5"], "pubchem":["pubchem5"],"lmid":["LMFA01010005"], "swissid" : ["swissid5"], "smiles":["smiles5"]},
                    "met6" : {"name":["lipid6"], "formula":["formula6"],"chebi":["00006"],"inchi":["1S/inchi6"],"hmdb":["hmdb6"], 
                                "kegg":["kegg6"], "pubchem":["pubchem6"],"lmid":["LMFA01010006"], "swissid" : ["swissid6"], "smiles":["smiles6"]}    
                }
     
    
    multimapping = mapping.multimapping(lipids, metabolites)
    
    
    metabolites_chebi_entities = mapping.convert_to_chebi_entities_once(metabolites)
    chebi_mapping = mapping.chebi_class_mapping(lipids, metabolites_chebi_entities)

    
    def test_multimapping_mapping(self):
        self.assertTrue("met1" in self.multimapping["lipid1"])
        
    def test_multimapping_type(self) :
        self.assertEqual(self.multimapping["lipid1"]["met1"]["mapping type"], "exact multimapping")

    def test_multimapping_ids(self) :
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["chebi"], ["00001"])
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["inchi"], ["1S/inchi1"])
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["lipidmaps"], ["LMFA01010001"])
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["kegg"], ["kegg1"])
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["pubchem"], ["pubchem1"])
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["hmdb"], ["hmdb1"])
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["exact smiles"], ["smiles1"])
        self.assertEqual(self.multimapping["lipid1"]["met1"]["identifiers"]["swisslipids"], ["swissid1"])
        
    
    def test_chebi_mapping(self):
        self.assertTrue("met3" in self.multimapping["lipid3"])

            
unittest.main()
