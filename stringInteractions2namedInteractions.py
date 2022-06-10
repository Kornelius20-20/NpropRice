"""

Author: Janith Weeraman

methods that takes a STRINGID formatted interaction file and then outputs an interaction file with names corresponding
to a given source or uniprot of the given source does not exist

"""

import gzip, os
import networkx as nx

aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
source_file = "txt/ppi.tsv"
output_file = "txt/named_ppi.tsv"
preferred_source = "Uniprot"
graphfile = 'test.gexf'


def create_aliasdict(aliasfile=aliasfile):
    """
    A method that takes a protein alias file as a gzip file and builds a dictionary with keys being STRINGIDs
    and values being dictionaries containing keys as sources and corresponding source-related names as values within
    the nested dict

    :param aliasfile: the protein alias file for an organism from STRING website
    :return: a dict with { STRINGID: { SOURCE : Name } } items
    """
    alias_key = {}  # dict to hold aliases

    with gzip.open(os.path.join(aliasfile), 'rt') as file:
        file.readline()
        line = file.readline()

        while line:
            # Split each line into it's items
            items = line.rstrip().split('\t')
            # Get the alias list for STRINGID of current line
            # Find the aliases for the stringID in the line from the dict and put into a temp variable
            tmp = alias_key.get(items[0], {})
            # get the source for the current line
            source = items[-1]
            # Add the new source and name into the alias list dict
            tmp[source] = items[1]
            # replace the aliasdict for the current STRINGID with the one with the current source added
            alias_key[items[0]] = tmp
            line = file.readline()

    return alias_key

def stringidconvert(proteinlist,aliasdict=create_aliasdict(aliasfile),source='BLAST_UniProt_ID'):

    outputproteins = []
    for prot in proteinlist:
        aliases = aliasdict.get(prot,None)

        if aliases is not None:
            outputproteins.append(aliases[source])
        else:
            outputproteins.append('')

    return outputproteins

def id2stringdict(aliasfile,source='BLAST_UniProt_ID'):
    """
    Function that will take in an aliasfile and convert it into an aliasdict dictionary. It will then create a new
    dictionary where the key is the name of a protein in the given source format andt he value is the STRING ID of the
    protein

    Parameters
    ----------
    aliasfile : the path to a file containing protein aliases
    source: the preferred source to have as a key. Default is BLAST_UniProt_ID

    Returns
    -------
    dict: dictionary of keys as source names and values as STRING IDs
    """

    alias_key = create_aliasdict(aliasfile)

    newdict = {}
    for key,item in alias_key.items():
        newkey = item[source]
        newdict[newkey] = key

    return newdict

if __name__ == "__main__":
    from cluster_drought_module_greedy import descendingdictkeys

    aliasdict = create_aliasdict(aliasfile)

    sources = {}
    numprots = len(aliasdict.keys())
    print(f"numprots: {numprots}")
    for item in aliasdict.values():
        for sc in item.keys():
            try:
                sources[sc] += 1
            except KeyError:
                sources[sc] = 1

    keys = descendingdictkeys(sources)
    print(keys)