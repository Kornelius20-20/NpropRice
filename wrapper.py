"""

Author: Janith Weeraman

Wrapper script to handle taking a graph and seedlist and running network propagation on it, then taking the output and
adding the GO terms and whether the node is a seed or not as metadata into the nodes. The resulting graphs will be
easier to analyze using a visualization software
Finally the script generates a csv file of the highest weighted proteins that are not seeds that resulted from the
network propagation

"""

import networkx as nx
import pandas as pd
import netprop
import approx_go
import os
from get_high_scores import get_non_seeds,descendingnodes

# From netprop
seedlist = r"txt/string_seeds.txt"
aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
dframe2 = "txt/uniprot_original.csv"
delim = ','
regen = False

# Load graph
maingraph = nx.read_gexf('graph.gexf')
# number of iterations to run to the random walk
iter = [25,30,35,40,45]
# following parameters should be given in list form
# weight to give to seeds
weight = [10]
# restart parameter
alpha = [0.1]

cutoffs = [25.0,30.0,35.0,40.0,45.0]

# Create directories to hold output files
if not os.path.exists("/outputs"):
    os.makedirs('/outputs')
if not os.path.exists("outputs/graphs"):
    os.makedirs('outputs/graphs')
if not os.path.exists("outputs/results"):
    os.makedirs('outputs/results')



dframe = pd.read_csv(dframe2,delimiter=delim)
frame = pd.read_csv(dframe2,delimiter=delim)

for cutoff in cutoffs:
    for i in range(len(weight)):
        for j in range(len(alpha)):
            for k in range(len(iter)):
                # Run network propagation with the given values
                graphfile = netprop.netprop(maingraph,seedlist,aliasfile,weight[i],alpha[j],iter[k],regen=regen,cutoff=cutoff)

                outgraph = graphfile[:-5] + "seedsAndGO.gexf"
                graph = nx.read_gexf(graphfile)


                graph = approx_go.assign_metadata(graph, frame, asterm=False,seedlist=seedlist)
                graph = approx_go.assign_best_go_id(graph)

                # Get nodes that aren't seeds and their weights
                nodelist, weightlist = get_non_seeds(graph, True)

                # sort them in ascending order
                nodes, weights = descendingnodes(nodelist, weightlist, True)

                entries = [dframe.loc[dframe['Entry'] == node] for node in nodes]
                if len(entries) > 0:
                    results = pd.concat(entries)
                    weights = [weights[nodes.index(i)] for i in results['Entry']]
                    results['weights'] = weights
                    results.to_csv(f'outputs/results/{graphfile[15:-5]}.csv')

                # nx.write_gexf(graph, outgraph)


import result_processing as rp

outputdir = "outputs/results"
outputprots = rp.concat_outputs(outputdir)
candset = set(outputprots['Entry'].tolist())
print(rp.get_best_csv_name(candset,outputdir))