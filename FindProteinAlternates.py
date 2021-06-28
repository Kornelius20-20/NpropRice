# -*- coding: utf-8 -*-
"""

Author: Janith Weeraman
Date: 28/06/2021
version: 1.0

A script to find all the different names and their sources for a given list of
proteins 


"""

from JCNUParser import aliasesFromFile

aliasfile = r"I:\Thesis\4530.protein.aliases.v11.0.txt.gz"

querylist = r"JCNUproteins.txt"
source = 'Ensembl_Uniprot'

# Create dataframe from alias file
frame = aliasesFromFile(aliasfile)

# Create list of proteins from a query file
with open(querylist,'r') as file:
    queryprot = [protein.strip() for protein in file.readlines()]

# For each protein in the qeury list
for protein in queryprot:
    # find its STRING ID and extract it
    stringID = frame.loc[frame['alias']==protein]['string_protein_id'].item()
    
    # Find all proteins with the same ID
    alternates = frame.loc[frame['string_protein_id']==stringID or frame['source'] == source]
    
    print(alternates[['alias','source']])