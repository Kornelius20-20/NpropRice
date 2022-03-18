import networkx as nx
import pandas as pd
import netprop
import approx_go

# From netprop
seedlist = r"txt/string_seeds.txt"
aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
regen = True

# Load graph
graph = nx.read_gexf('graph.gexf')
# number of iterations to run to the random walk
iter = [50]
# following parameters should be given in list form
# weight to give to seeds
weight = [100]
# restart parameter
alpha = [0.5,1.2,10,0.001]
# Run network propagation with the given values
graphfiles = netprop.netprop(graph,seedlist,aliasfile,weight,alpha,iter,regen=True)
print(graphfiles)
# graphfiles = ['outputs/p-100-0.5-50.gexf', 'outputs/p-100-1.2-50.gexf', 'outputs/p-100-10-50.gexf',
#               'outputs/p-100-0.001-50.gexf']

# From approxgo.py
# graphfiles = ["p-1000-0.5-100.gexf"]
dframe2 = "txt/uniprot_original.csv"

# from get_high_scores.py
from get_high_scores import get_non_seeds,descendingnodes
import pandas as pd
dframe = pd.read_csv("txt/uniprot_original.csv")
frame = pd.read_csv(dframe2)
for graphfile in graphfiles:
    outgraph = graphfile[:-5] + "seedsAndGO.gexf"
    graph = nx.read_gexf(graphfile)


    graph = approx_go.assign_metadata(graph, frame, asterm=False)
    graph = approx_go.assign_best_go_id(graph)


    # graphfiletoread = graphfile[:-5] + "seedsAndGO.gexf"
    print(graphfile)
    # graph = nx.read_gexf(graphfiletoread)

    # Get nodes that aren't seeds and their weights
    nodelist, weightlist = get_non_seeds(graph, True)

    # sort them in ascending order
    nodes, weights = descendingnodes(nodelist, weightlist, True)

    entries = [dframe.loc[dframe['Entry'] == node] for node in nodes]
    if len(entries) > 0:
        results = pd.concat(entries)
        weights = [weights[nodes.index(i)] for i in results['Entry']]
        results['weights'] = weights
        results.to_csv(f'outputs/results/{graphfile[8:]}.csv')

    nx.write_gexf(graph, outgraph)

from playsound import playsound
playsound("test.mp3")