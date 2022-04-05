import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import approx_go
import pandas as pd


dframe2 = "txt/uniprot_original.tsv"
graph = nx.read_gexf('outputs/graphs/p-10-0.1-50seedsAndGO.gexf')
frame = pd.read_csv(dframe2,delimiter='\t')

graph = approx_go.assign_metadata(graph, frame, asterm=False)

results = greedy_modularity_communities(graph)

for i in range(len(results)):
       for node in results[i]:

              try:
                     print(graph.nodes[node]['GO'])
              except KeyError:
                     continue



# nx.write_gexf(graph,'outputs/graphs/p-10-0.1-50seedsAndGO-w-community.gexf')