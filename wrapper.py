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
from newnetprop import netprop
import os

import stringInteractions2namedInteractions
from get_high_scores import get_non_seeds,descendingnodes

# From netprop
seedlist = r"txt/string_seeds.txt"
aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
dframe2 = "txt/uniprot_original.csv"
delim = ','
outputdir = "outputs/graphs"
regen = False

# Load graph
maingraph = nx.read_gexf('graph_test.gexf')
# number of iterations to run to the random walk
iter = [5]
# following parameters should be given in list form
# weight to give to seeds
weight = [100]
# restart parameter
alpha = [0.1]

# cutoffs = [10.0]

# Create directories to hold output files
if not os.path.exists("/outputs"):
    os.makedirs('/outputs')
if not os.path.exists("outputs/graphs"):
    os.makedirs('outputs/graphs')
if not os.path.exists("outputs/results"):
    os.makedirs('outputs/results')



dframe = pd.read_csv(dframe2,delimiter=delim)
frame = pd.read_csv(dframe2,delimiter=delim)

def generate_graphs(weight,alpha,iter):
    for i in range(len(weight)):
        for j in range(len(alpha)):
            for k in range(len(iter)):
                # Run network propagation with the given values
                graphfile = netprop.netprop(maingraph, seedlist, weight[i], alpha[j], iterevery=5,iterend=10, regen=regen)

def cutoff_graph(graphname,attr,cutoffs):

    graph = nx.read_gexf(os.path.join("outputs",graphname))

    for cut in cutoffs:
        nodes = []

        for node in graph.nodes:
            if graph.nodes[node][attr] > cut: nodes.append(node)

        newgraph = nx.subgraph(graph,nodes)

        nx.write_gexf(newgraph,os.path.join(outputdir,f"{graphname[:-5]}_{cut}.gexf"))

def get_top(graph,attr,howmany=100,returnscores=False):
    from get_high_scores import descendingnodes

    keys, values = zip(*nx.get_node_attributes(graph,attr).items())

    topnamesandscores = descendingnodes(keys,values,returnboth=returnscores)

    if returnscores:
        topnamesandscores[0] = topnamesandscores[0][:howmany]
        topnamesandscores[1] = topnamesandscores[1][:howmany]

        return (topnamesandscores[0], topnamesandscores[1])
    else:
        topnamesandscores = topnamesandscores[:howmany]

        return topnamesandscores


    return graph

def add_seeds(graph,seedlist):
    # Mark nodes
    for node in seedlist:
        try:
            graph.nodes[node]['isSeed'] = True
        except KeyError:
            None

    return graph


if __name__=="__main__":
    from stringInteractions2namedInteractions import stringidconvert,create_aliasdict

    aliasdict = create_aliasdict(stringInteractions2namedInteractions.aliasfile)

    manualseeds = "txt/string_seeds.txt"
    with open(manualseeds, 'r') as file:
        manualseeds = [line.strip() for line in file.readlines()]


    # generate_graphs(weight,alpha,iter)

    iterations = [50,100,250,500]

    for i in iterations:

        for _,_,files in os.walk("outputs"):

            for file in files:
                graph = nx.read_gexf(os.path.join("outputs",file))

                predicts = get_top(graph,'propagated_weight',howmany=i)
                predicts.extend(manualseeds)

                graph = nx.subgraph(graph,predicts)
                graph = add_seeds(graph,manualseeds)

                nx.write_gexf(graph,os.path.join("outputs/graphs",f"{file[:-5]}_{i}.gexf"))

            break
