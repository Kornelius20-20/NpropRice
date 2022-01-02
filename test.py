import networkx as nx
import pandas as pd


graphfile = "graph.gexf"
dframe = "txt/processed_uniprot.csv"

# graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)


import re

for gofile in frame["Gene ontology (biological process)"].to_list():
    if not type(gofile) is float:
        res = re.findall("GO:\d\d\d\d\d\d\d",gofile)
        print(res)