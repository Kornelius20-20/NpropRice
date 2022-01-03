import networkx as nx
import pandas as pd
import re

graphfile = "drought_module.gexf"
dframe = "txt/processed_uniprot.csv"

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

# Add gene ontology IDs and reviewed status as attributes
nodes = {}
for prot in graph.nodes:
        status = frame[frame['Entry'] == prot]
        if status.empty:
            None
        else:
            attributes = {'status':status['Status'].to_string()[-8:]}

            nodes[prot] = attributes
nx.set_node_attributes(graph, nodes)


# write graph as gexf file
nx.write_gexf(graph, "drought_module2.gexf")