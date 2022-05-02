def rescale(valuedict,scalemax=256,mincutoff=0):
    """
    Method to rescale a dictionary with numeric values

    :param valuedict: the dictionary with keys and numeric values
    :param scalemax: the max value to scale to
    :param mincutoff: the value below which inputs will be rounded to 0
    :return:
    """

    scalemax = scalemax -1
    smol = min(valuedict.values())
    big = max(valuedict.values())
    for key in valuedict.keys():
        newval = int(valuedict[key]/(big-smol) * scalemax)
        valuedict[key] = newval if newval >= mincutoff else 0

"""

Using the GO REST API to get the related GO terms for a particular term

"""

def get_go_json(goterm):
    import requests, sys, json
    requestURL = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{goterm}/descendants?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"

    r = requests.get(requestURL, headers={"Accept": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.json()

    with open('response.json', 'w') as file:
        json.dump(responseBody, file)

def change_nodes_of_gexf(graphfile,alias_key, preferred_source,backup_source= "Uniprot"):
    import networkx as nx

    graph = nx.read_gexf(graphfile)

    for node in graph.nodes:
        graph.nodes[node]['label'] = alias_key[node].get(preferred_source, alias_key[node][backup_source])

    nx.write_gexf(graph,'test2.gexf')


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

# code that formats the clusters to be in a tsv format with a number of columns equal to the
                     # longest cluster
                     # longestlist = max([len(item) for item in proteins])
                     #
                     # with open(os.path.join(outputdir, f'clusterfile{file[:-5]}.tsv'), 'w') as multlst:
                     #        multlst.write('\t'.join(titles))
                     #        multlst.write('\n')
                     #
                     #        for i in range(longestlist):
                     #               line = ''
                     #               for clust in proteins:
                     #                      try:
                     #                             line += clust[i] + '\t'
                     #                      except IndexError:
                     #                             line += '\t'
                     #               line += '\n'
                     #               multlst.writelines(line)

if __name__ == "__main__":
    goterm = "GO:0009819" # drought response
    filename = "response.json"
    get_go_json(goterm)