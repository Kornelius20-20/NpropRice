import os
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

def greedyclustergraph(graph, frame, aliasfile,id='BLAST_UniProt_ID',asterm=False):

       # graph = approx_go.assign_metadata(graph, frame, asterm=asterm)

       results = greedy_modularity_communities(graph)

       for i in range(len(results)):
              # add cluster label
              for node in results[i]:
                     graph.nodes[node]['cluster'] = i + 1


       from stringInteractions2namedInteractions import stringidconvert,create_aliasdict

       aliasdict = create_aliasdict(aliasfile)

       proteins = [stringidconvert(results[i],aliasdict,id) for i in range(len(results))]

       return proteins


aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
outputdir = "outputs/results"
dframe2 = "txt/uniprot_original.csv"
frame = pd.read_csv(dframe2,delimiter=',')

for _,_,files in os.walk('outputs/graphs'):
       for file in files:
              if "seedsAndGO" not in file:
                     graph = nx.read_gexf(os.path.join('outputs/graphs',file))

                     proteins = greedyclustergraph(graph,frame,aliasfile)

                     titles = [f"cluster{i + 1}" for i in range(len(proteins))]
                     longestlist = max([len(item) for item in proteins])

                     with open(os.path.join(outputdir, f'clusterfile{file[:-5]}.tsv'), 'w') as multlst:
                            multlst.write('\t'.join(titles))
                            multlst.write('\n')

                            for i in range(longestlist):
                                   line = ''
                                   for clust in proteins:
                                          try:
                                                 line += clust[i] + '\t'
                                          except IndexError:
                                                 line += '\t'
                                   line += '\n'
                                   multlst.writelines(line)

                     nx.write_gexf(graph,os.path.join('outputs/graphs',file))
