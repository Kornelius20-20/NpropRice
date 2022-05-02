import os
import networkx as nx
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


# calculating partition coefficient
def partition_coefficient(graph, cluster_attr='cluster'):
       # Get total number of clusters
       totclust = max(nx.get_node_attributes(graph, cluster_attr).values())

       # for each node
       for node in graph.nodes:
              module_part = 0
              # For each module
              for i in range(1, totclust + 1):
                     # Get its neightbors and find out how much of them are belong to the current cluster of the iteration
                     nb_counter = 0
                     for neighbor in nx.neighbors(graph, node):
                            if graph.nodes[neighbor][cluster_attr] == i:
                                   nb_counter += 1

                     # sum of partitions of each module
                     try:
                            module_part += (nb_counter / nx.degree(graph, node)) ** 2
                     except ZeroDivisionError:
                            module_part = module_part
              part_coef = 1 - module_part

              graph.nodes[node][f'{cluster_attr}_PC'] = part_coef

       return graph


def get_best_scoring_nodes(graph, attr, cutoff):
       value_attrs = nx.get_node_attributes(graph, attr)

       # resccle the attribute values
       maxval = max(value_attrs.values())
       minval = min(value_attrs.values())
       for key in value_attrs.keys():
              value_attrs[key] = (value_attrs[key] - minval) / (maxval - minval) * 100

       # Only get values larger than the cutoff
       outputnodes = []
       for key, value in value_attrs.items():
              if value >= cutoff: outputnodes.append([value, key])

       outputnodes.sort(reverse=True)

       return [i[1] for i in outputnodes]


def add_clustering_as_attr(graph, clusters, attrname):
       clust = 1
       for cluster in clusters:
              for node in cluster:
                     graph.nodes[node][attrname] = clust
              clust += 1

       return graph


aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
outputdir = "outputs/results"
dframe2 = "txt/uniprot_original.csv"
frame = pd.read_csv(dframe2,delimiter=',')

for _,_,files in os.walk('outputs/graphs'):
       for file in files:
              if "seedsAndGO" not in file:
                     graph = nx.read_gexf(os.path.join('outputs/graphs',file))

                     # Do greedy clustering
                     attr = 'greedy_clusters'
                     proteins = nx.community.greedy_modularity_communities(graph)

                     graph = add_clustering_as_attr(graph,proteins,attr)
                     graph = partition_coefficient(graph,attr)

                     # Do label propagation
                     attr = 'label_propagation'
                     proteins = nx.community.label_propagation_communities(graph)

                     graph = add_clustering_as_attr(graph, proteins, attr)
                     graph = partition_coefficient(graph, attr)


              nx.write_gexf(graph,os.path.join('outputs/graphs',file))
