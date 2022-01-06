import networkx as nx
import pandas as pd
import re
from collections import Counter


graphfile = "graph.gexf"
dframe = "txt/processed_uniprot.csv"

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

pd.set_option('display.max_colwidth', 1000)


# Add GO terms from the GO columns into a 'GO' attribute as a set
for node in graph.nodes:
    row = frame[frame['Entry'] == node]

    if not row.empty:
        # For every node who is present in the uniprot data, get GO IDs
        biogo = row["Gene ontology (biological process)"].to_string() # + row["Gene ontology (GO)"].to_string()
        # biogo += row["Gene ontology (molecular function)"].to_string() + row["Gene ontology IDs"].to_string()
        # Extract just the GO terms and form them into a non-redundant set using "GO:\d\d\d\d\d\d\d" regex
        firstres = re.findall("    .*? \[GO", biogo)
        res = re.findall("; .*? \[GO", biogo)
        res = [item[2:-4] for item in res]
        firstres = [item[4:-4] for item in firstres]

        goset = set(res + firstres)



        # Add to the 'GO' attribute
        if len(goset) > 0:
            graph.nodes[node]['GO'] = goset
            graph.nodes[node]['isSeed'] = True


# Iterate over all nodes, getting nodes and their attributes
for node,attrs in graph.nodes(True):
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


nx.write_gexf(graph,"testgraph.gexf")
