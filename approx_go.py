"""

Script to take the formed PPI network and assign metadata values to the nodes such as GO IDs, protein names, whether the
given node is a present in the seed list etc.

"""

import networkx as nx
import pandas as pd
from collections import Counter


graphfiles = ["p -100-0.5.gexf","p -10-2.gexf"]
dframe = "txt/processed_uniprot.csv"



def assign_metadata(graph,infoframe,asterm=True):
    import re

    # Increase column width of pd output
    pd.set_option('display.max_colwidth', 1000)
    # Add GO terms from the GO columns into a 'GO' attribute as a set
    for node in graph.nodes:
        row = infoframe[infoframe['Entry'] == graph.nodes[node]['label']]

        if not row.empty:
            # For every node who is present in the uniprot data, get GO IDs
            biogo = row["Gene ontology (biological process)"].to_string() # + row["Gene ontology (GO)"].to_string()
            # biogo += row["Gene ontology (molecular function)"].to_string() + row["Gene ontology IDs"].to_string()

            if asterm:
                # Extract just the GO terms and form them into a non-redundant set using "GO:\d\d\d\d\d\d\d" regex
                firstres = re.findall("GO:\d\d\d\d\d\d\d", biogo)
                goset = set(firstres)
            else:
                # Extract the GO terms and add their associated names only to the list
                firstres = re.findall("    .*? \[GO", biogo)
                res = re.findall("; .*? \[GO", biogo)
                res = [item[2:-4] for item in res]
                firstres = [item[4:-4] for item in firstres]
                goset = set(res + firstres)



            # Add to the 'GO' attribute
            if len(goset) > 0:
                graph.nodes[node]['GO'] = goset
                graph.nodes[node]['isSeed'] = True

                # from pullingGOviaREST.py. Adds protein name in place of node label
                name = row['Protein names'].to_string()

                result = re.findall("\(.*?\)", name)

                # If there's a result, replace label with last result
                if len(result) > 0: graph.nodes[node]['label'] = result[-1][1:-1]

    return graph

def assign_best_go_id(graph):
    # Iterate over all nodes, getting nodes and their attributes
    for node, attrs in graph.nodes(True):
        # For nodes with GO terms in attributes
        if 'GO' in attrs:
            # Get the GO terms
            nodego = attrs['GO']
            combinedgos = []

            # For each neighbor of node
            for neighbor in graph.neighbors(node):
                # Try to get the neighbors GO terms if it exists
                try:
                    neighborgo = graph.nodes[neighbor]['GO']
                except KeyError:
                    continue

                # Find GO terms that are in common with neighbor and add to combinedgos list
                intersected = nodego.intersection(neighborgo)
                combinedgos += list(intersected)

            # Count the occurences of each GO term in the list and assign 'GO' value as highest value in list
            golabels = Counter(combinedgos)
            mostcomlabel = golabels.most_common(1)
            if len(mostcomlabel) > 0:
                graph.nodes[node]['GO'] = mostcomlabel[0][0]

            else:
                graph.nodes[node]['GO'] = list(graph.nodes[node]['GO'])[0]

    return graph


if __name__ == "__main__":
    for graphfile in graphfiles:
        outgraph = graphfile[:-5] + "seedsAndGO.gexf"
        graph = nx.read_gexf(graphfile)
        frame = pd.read_csv(dframe)

        graph = assign_metadata(graph,frame,asterm=False)
        graph = assign_best_go_id(graph)
        nx.write_gexf(graph,outgraph)
