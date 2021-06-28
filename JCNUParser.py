"""

Author: Janith Weeraman
Date created: 26/06/2021
Date modified: 28/06/2021
Version: 1.2

Script to get the genes from drought gene databse in http://ccbb.jnu.ac.in
and convert the gene names into a format that can be procssed with string 
database files.

In order to get for a different organism, enter the relavant organism genus 
name in the 'organism' variable

NOTE:Orysa sativa -> oryza (only this organism needs to have a non-capitalized
                            genus name)
     Arabidopsis thaliana -> Arabidopsis


    
"""

import requests
from bs4 import BeautifulSoup
import pandas as pd

# Organism in the database for which you want the stress genes
organism = "oryza" 
# Page of JNU stress gene database
database_site = f"http://ccbb.jnu.ac.in/stressgenes/{organism}.html" 
# Location of protein alias file for your organism
aliasfile = r"I:\Thesis\4530.protein.aliases.v11.0.txt.gz" 


def aliasesFromFile(aliasfile):
    """
    Function to return a set of aliases from a gz file of protein aliases for 
    an organism
        

    Parameters
    ----------
    aliasfile : String
        A link to the alias file given as a String.

    Returns
    -------
    frame : pd.DataFrame
        A dataframe built from the alias file with colums ['string_protein_id',
                                                           'alias', 'source'].

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
    
    # Decode the document
    data = [line.decode('utf-8').strip() for line in data]


    # Get only titles from the header line
    titles = data[0].split("##")[1:-1]
    titles = [title.strip() for title in titles]
    
    
    # Split data lines into lists of values separated by tab delimitter
    data = [line.split('\t') for line in data[1:]]
    
    
    frame = pd.DataFrame(data,columns=titles)
    
    return frame


def proteinsInAlias(proteins,aliasfile):
    """
    Function to return a set of genes that are common between a list of give
    n proteins and a list of proteins generated from a protein alias file from
    STRING db

    Parameters
    ----------
    proteins : list
        iterable of proteins as strings.
    aliasfile : String
        link to path of alias file given as a string.

    Returns
    -------
    TYPE set
        A set of proteins that are in common between database and alias file.

    """
    
    # Create dataframe from alias file
    aliases = aliasesFromFile(aliasfile)
    
    proteins = [protein.upper() for protein in proteins]
    
    # Find common proteins between two sets    
    return set(proteins).intersection(set(aliases['alias'])) 


def main():
    """
    
    If script is run as main then the function will get the stress genes
    associated with Oryza sativa from the database and find which ones
    are in common with a given STRING alias file

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Get stress gene databse
    page = requests.get(database_site)

    # Throw an error if the page get fails
    if page.status_code == 404: 
        raise Exception("""Couldn't find gene database for organism.
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

