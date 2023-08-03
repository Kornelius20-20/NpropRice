# -*- coding: utf-8 -*-
"""

Author: Janith Weeraman
Date: 28/06/2021
version: 1.0

A script to find all the different names and their sources for a given list of
proteins 


"""

from oldfiles.JCNUParser import aliasesFromFile

aliasfile = r"gz/4530.protein.aliases.v11.0.txt.gz"

querylist = r"JCNUproteins.txt"
source = 'Ensembl_UniProt'

def proteinNameConverter(queryprot,source=source,aliasfile=aliasfile):
    """
    A function to take a list of proteins and find all alternative names for 
    the proteins that exist in the Uniprot database, using a STRING alias file

    Parameters
    ----------
    queryprot : list
        A list of protein names as strings.
    source : string
        What protein output type do you want. default: Ensembl_UniProt.
    aliasfile : String, optional
        path to a protein alias file. The default is a hardcoded test path.

    Returns
    -------
    outputdict : dict
        A dictionary with { query_name:[list of ids] } pairs.

    """

    # Create dataframe from alias file
    frame = aliasesFromFile(aliasfile)
        
    outputdict = {}
    
    # For each protein in the qeury list
    for protein in queryprot:
        # find its STRING ID and extract it
        stringID = frame.loc[frame['alias']==protein]['string_protein_id'].item()
        
        # Find Uniprot ID of proteins with the same ID
        alternates = frame.loc[(frame['string_protein_id']==stringID) & (frame['source'] == source)]
        
        # Create a dict for holding name: [uniprot_IDs] pairs
        outputdict[protein] = []
        
        for item in alternates['alias']:
            outputdict[protein].append(item)
            
    
    return outputdict


if __name__=="__main__":
    
    
    # Create list of proteins from a query file
    with open(querylist,'r') as file:
        queryprot = [protein.strip() for protein in file.readlines()]
    
    print(proteinNameConverter(queryprot, source))
