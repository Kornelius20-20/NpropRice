import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import approx_go
import pandas as pd

def descendingdictkeys(dic,desc=True):
       """
       Method that takes in a dictionary with keys and values that are integers and returns the keys in descending order
       of the values

       Parameters
       ----------
       dic: input dictionary

       Returns
       -------
       list of keys sorted based on order specified

       """

       keys, values = zip(*dic.items())
       keys = list(keys)
       values = list(values)

       output = []
       while len(values) > 0:
              bestind = values.index(max(values))
              output.append(keys[bestind])
              keys.pop(bestind)
              values.pop(bestind)

       if desc:
              return output
       else: return output.reverse()


dframe2 = "txt/uniprot_original.tsv"
graph = nx.read_gexf('outputs/graphs/p-10-0.1-50seedsAndGO.gexf')
frame = pd.read_csv(dframe2,delimiter='\t')

graph = approx_go.assign_metadata(graph, frame, asterm=False)

results = greedy_modularity_communities(graph)

for i in range(len(results)):

       # for each node add the sum of the degrees of neighbor nodes
       for node in results[i]:
              neighbors = nx.neighbors(graph,node)
              degreesum = 0
              for nb in neighbors:
                     degreesum += nx.degree(graph,nb)
              graph.nodes[node]['neighbordegrees'] = degreesum


       # add cluster label
       for node in results[i]:
              graph.nodes[node]['cluster'] = i + 1

       gos_in_clust = {}

       # For nodes with GO terms, score their terms in the cluster GO list
       for node in results[i]:
              try:
                     assignedgo = graph.nodes[node]['GO']
                     if type(assignedgo) is str:
                            score = gos_in_clust.get(assignedgo,0)
                            gos_in_clust[assignedgo] = score + 1
                     else:
                            for go in assignedgo:
                                   score = gos_in_clust.get(go, 0)
                                   gos_in_clust[go] = score + 1

              except KeyError:
                     continue

       orderedgos = descendingdictkeys(gos_in_clust)

       for node in results[i]:
              try:
                     del graph.nodes[node]['GO']
                     for i in range(0,5):
                            goterm = f'GO{i+1}'
                            graph.nodes[node][goterm] = orderedgos[i]

              except KeyError:
                     None


nx.write_gexf(graph,'outputs/graphs/p-10-0.1-50seedsAndGO-w-community.gexf')