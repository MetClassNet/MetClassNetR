# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 17:07:34 2018

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
#coding: "utf8"

import libchebipy as libchebi
"""
Warning libchebi needs str to create a chebi entity, but will return chebi id as int
This is why in most functions we change the type of the return to str so that we can easily use it in further functions
"""

def normalize_chebi_id(chebi):
    """
    Removes the ChEBI: in the chebi id strings, will be used so that every ChEBI is only a string of number for consistency
    Args:
        arg1 (str): chebi identifier
    Returns:
        str: chebi identifier, but only the part made of numbers
        
    """
    if chebi != "" and chebi != None:
        if "chebi:" in chebi.lower():
            chebi = chebi[6:]
            return chebi
        else:
            return chebi
    
    return chebi


def get_true_id(chebi_entity):
    """
    Get the main id of a given chebi id (molecules can have multiple chebi id, one of them is the main one)
    Important to ensure we compare the same identifiers. 
    For example a molecule can be described in different files with different identifiers that are in fact linked to the same amin identifier
    
    Args: 
        arg1:ChEBI entity
    
    Returns: 
        str: Main ChEBI ID of the entity argument
    """
    
    true_id = chebi_entity.get_parent_id()
    if true_id != None:
        return true_id
    else:
        return chebi_entity.get_id()
    

def get_tautomer(chebi_entity):
    """
    Check if the molecule has a conjugate acid or base, or a tautomer using chebi identifiers,
    and returns it if there is one
    
    Args:
        arg1: ChEBI entity (obtained with libchebi)
    Returns: 
        ChEBI entity: ChEBI entity of tautomer if there is one, else None
    
    """
    tautomers_id = []
    outgoings = chebi_entity.get_outgoings()
    for relation in outgoings:
        if relation.get_type() == "is_tautomer_of":
            tautomer_id = relation.get_target_chebi_id()
            if tautomer_id != "":
                true_tautomer_id = get_true_id(libchebi.ChebiEntity(tautomer_id))
                tautomers_id.append(normalize_chebi_id(true_tautomer_id))
            #return true_tautomer_id
    
    return tautomers_id


def get_conjugate_acid_base(chebi_entity):
    """
    Check if the molecule has a conjugate acid or base, or a tautomer using chebi identifiers,
    and returns it if there is one
    
    Args: 
        arg1: ChEBI entity
    Returns: 
        ChEBI entity: of conjugate base or acid, if there is one, else None
    """
    conjugates_id = []
    outgoings = chebi_entity.get_outgoings()
    for relation in outgoings:
        if relation.get_type() == "is_conjugate_acid_of" or relation.get_type() == "is_conjugate_base_of":
            conjugate_id = relation.get_target_chebi_id()
            if conjugate_id != "":
                true_conjugate_id = get_true_id(libchebi.ChebiEntity(conjugate_id))
                conjugates_id.append(normalize_chebi_id(true_conjugate_id))

    return conjugates_id

def get_parents(chebi_entity):
    """
    Use the methods of libchebi to get the parents of a ChEBI entity (molecule)
    
    Args: 
        arg1: ChEBI entity of molecule
    Returns: 
        list: of ChEBI entities that are parents to the query
    """
    outgoings = chebi_entity.get_outgoings()
    parents = []
    for relation in outgoings:
        if relation.get_type() == "is_a":
            parent_id = relation.get_target_chebi_id()
            parent = libchebi.ChebiEntity(parent_id)
            parents.append(parent)
    
    return parents


def get_children(chebi_entity):
    """
    Use the methods of libchebi to get the children of a ChEBI entity (molecule)

    Args:
        arg1: ChEBI entity of molecule
    Returns:
        list: of ChEBI entities that are children to the query
    """

    incomings = chebi_entity.get_incomings()
    children = []
    for relation in incomings:
        if relation.get_type() == "is_a":
            child_id = relation.get_target_chebi_id()
            child = libchebi.ChebiEntity(child_id)
            children.append(child)

    return children


def create_parent_graph(chebi_entity, graph):
    """
    Recursively recreate the ontology graph containing the PARENTS of the molecule using chebilib
    
    Args: Chebi id of molecule
          arg1 : empty graph that will be built recursively
    
    Returns:
        node : Parent of the molecule used for this recursion, will be used as base in next recursion
        graph: Graph of the ontology (dictionary of dictionaries, each node as key has its linked nodes in values)
    """
    chebi_id = normalize_chebi_id(get_true_id(chebi_entity))

    parents = get_parents(chebi_entity)
    parents_id = []
    for parent in parents:
        parent_id = normalize_chebi_id(get_true_id(parent))
        parents_id.append(parent_id)
    
    node = {chebi_id : parents_id}
    graph.update(node)
    
    for parent_id in parents_id:
        parent_entity = libchebi.ChebiEntity(parent_id)
        create_parent_graph(parent_entity,graph)

    return node,graph


def create_children_graph(chebi_entity, graph):
    """
    Recursively recreate the ontology graph containing the CHILDREN of the molecule using chebilib

    Args: Chebi id of molecule
          arg1 : empty graph that will be built recursively

    Returns:
        node : Parent of the molecule used for this recursion, will be used as base in next recursion
        graph: Graph of the ontology (dictionary of dictionaries, each node as key has its linked nodes in values)
    """
    chebi_id = normalize_chebi_id(get_true_id(chebi_entity))

    children = get_children(chebi_entity)
    children_id = []
    for child in children:
        child_id = normalize_chebi_id(get_true_id(child))
        children_id.append(child_id)

    node = {chebi_id: children_id}
    graph.update(node)

    for child_id in children_id:
        child_entity = libchebi.ChebiEntity(child_id)
        create_children_graph(child_entity, graph)

    return node, graph



def find_shortest_path(graph, start, end, path=[]):
    """
    Recusrsively finds the shortest path in the graph created with create_graph between two chebi IDs
    Args:
        arg1: Graph to look into
        arg2 (str): Start chebi
        arg3 (str): End chebi
        arg4 (list): Path empty path that will used  in the recursion to store the shortest path
    Returns:
        list: Shortest path, list of chebi identifiers from start to end
    """
    start = normalize_chebi_id(start)
    end = normalize_chebi_id(end)
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = find_shortest_path(graph, node, end, path)
        if newpath != None:
            if shortest == None or len(newpath) < len(shortest):
                shortest = newpath
    return shortest



