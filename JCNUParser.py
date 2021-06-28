"""

Author: Janith Weeraman
Date: 26/06/2021

Script to get the genes from drought gene databse in http://ccbb.jnu.ac.in
and convert the gene names into a format that can be procssed with string database files.

In order to get for a different organism, enter the relavant organism genus name in the 
'organism' variable

NOTE:Orysa sativa -> oryza (only this organism needs to have a non-capitalized genus name)
     Arabidopsis thaliana -> Arabidopsis


    
"""

import requests
from bs4 import BeautifulSoup
import pandas as pd


organism = "oryza" # Organism in the database for which you want the stress genes
database_site = f"http://ccbb.jnu.ac.in/stressgenes/{organism}.html" # Page of JNU stress gene database
aliasfile = r"I:\Thesis\4530.protein.aliases.v11.0.txt.gz" # Location of protein alias file for your organism


def aliasesFromFile(aliasfile):
    """
    Function to return a set of aliases from a gz file of protein aliases for an organism
    
    aliasfile : String link
    """
    
    import gzip
    
    # load protein alias data for organism
    try:
        with gzip.open(aliasfile,'r') as file:
            data = file.readlines()
    except FileNotFoundError:
        aliasfile = input("Path to protein alias file: ")
        with gzip.open(aliasfile,'r') as file:
            data = file.readlines()
        
    aliases = set() # New set for holding only protein aliases
    for datum in data[1:]:
        # Decode the line, process it and split it using tab delimitters
        line = datum.decode('utf-8').strip().split('\t')
        aliases.add(line[1])
    
    return aliases


def proteinsInAlias(proteins,aliasfile):
    """
    
    Function to return a set of genes that are common between a list of givn proteins and
    a list of proteins generated from a protein alias file from STRING db
    
    proteins: iterable (list,set etc)
    aliasfile: String link
    """
    
    aliases = aliasesFromFile(aliasfile)
        
    return set(proteins).intersection(aliases) # Find common proteins between two sets


def main():
    # Get stress gene databse
    page = requests.get(database_site)

    # Throw an error if the page get fails
    if page.status_code == 404: 
        raise Exception(f"""Couldn't find the stress gene database for {organism}.
        Try changing the organism name or check viable organisms in 
        http://ccbb.jnu.ac.in/stressgenes/""")
        input()
    
    
    # Convert page to soup object and create a dataframe of the gene table
    soup = BeautifulSoup(page.content,'html.parser')
    table = soup.find_all('table')[2:] # Get data table without caption

    frame = pd.read_html(str(table))[0]
    
    genes = proteinsInAlias(frame['Gene Name'],aliasfile)
    
    print('Genes common between JNU stress database and STRING db aliases: \n')
    for gene in genes:
        print(gene)


if __name__ == '__main__':
    main()
    input()

