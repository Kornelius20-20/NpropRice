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


def create_aliasdict(aliasfile):
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


def id_translate(source_file, output_file, alias_key, preferred_source,header=False):
    """
    A method that takes in a source file of interactions and converts it into an output tsv containing the protein
    interactions with names as specified by preferred_source. If the preferred source does not exist then Uniprot is
    used as a fallback source

    :param source_file:
    :param output_file:
    :param alias_key:
    :return:
    """
    backup_source = "Uniprot"
    out = open(os.path.join(output_file), 'w')

    with open(os.path.join(source_file), 'r') as file:
        line = file.readline()
        if header: line = file.readline()

        while line:
            items = line.rstrip().split('\t')
            # For each STRING_ID, get the protein name from the preferred source in the alias file
            # If that is missing then go for the backup source
            left_prot = alias_key[items[0]].get(preferred_source, alias_key[items[0]][backup_source])
            right_prot = alias_key[items[1]].get(preferred_source, alias_key[items[1]][backup_source])
            out.writelines('\t'.join([left_prot, right_prot, items[-1], '\n']))
            line = file.readline()

    out.close()

def change_nodes_of_gexf(graphfile,alias_key, preferred_source,backup_source= "Uniprot"):
    graph = nx.read_gexf(graphfile)

    for node in graph.nodes:
        graph.nodes[node]['label'] = alias_key[node].get(preferred_source, alias_key[node][backup_source])

    nx.write_gexf(graph,'test2.gexf')

if __name__ == "__main__":
    alias_key = create_aliasdict(aliasfile)
    # id_translate(source_file, output_file, alias_key, preferred_source,header=True)
    change_nodes_of_gexf(graphfile,alias_key,'Uniprot')
