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
import netprop_copy as netprop
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
maingraph = nx.read_gexf('graph.gexf')
# number of iterations to run to the random walk
iter = [40]
# following parameters should be given in list form
# weight to give to seeds
weight = [10]
# restart parameter
alpha = [0.1,1.1,2]

cutoffs = [float(i) for i in range(50,76,5)]

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
                graphfile = netprop.netprop(maingraph, seedlist, aliasfile, weight[i], alpha[j], iterevery=10, regen=regen)

def cutoff_graph(graphname,attr,cutoffs):

    graph = nx.read_gexf(os.path.join("outputs",graphname))

    for cut in cutoffs:
        nodes = []

        for node in graph.nodes:
            if graph.nodes[node][attr] > cut: nodes.append(node)

        newgraph = nx.subgraph(graph,nodes)

        nx.write_gexf(newgraph,os.path.join(outputdir,f"{graphname[:-5]}_{cut}.gexf"))

def get_top_100(graphname,attr,howmany=100):
    from get_high_scores import descendingnodes

    graph = nx.read_gexf(os.path.join("outputs", graphname))

    keys, values = zip(*nx.get_node_attributes(graph,attr).items())

    top = descendingnodes(keys,values)[:howmany]

    return top

# generate_graphs(weight,alpha,iter)

for _,_,files in os.walk("outputs"):

    for file in files:
        cutoff_graph(file,'weight',cutoffs)

    break

def get_drought_module(graphfile,nodes):

    nodes = set(nodes)
    graph = nx.Graph()

    with open(graphfile,'r') as file:
        line = file.readline()
        while True:
            line = file.readline()
            if line == '': break

            line = line.strip().split('\t')
            if line[0] in nodes and line[1] in nodes:
                graph.add_edge(line[0],line[1],weight=line[2])


    return graph


from stringInteractions2namedInteractions import stringidconvert,create_aliasdict

aliasdict = create_aliasdict(stringInteractions2namedInteractions.aliasfile)

manualseeds = "txt/string_seeds.txt"
with open(manualseeds, 'r') as file:
    manualseeds = [line.strip() for line in file.readlines()]


for _,_,files in os.walk("outputs/graphs"):

    for file in files:

        graph = nx.read_gexf(os.path.join("outputs/graphs",file))

        # Get predicted drought proteins
        preds = list(graph.nodes)
        # add seeds
        preds.extend(manualseeds)

        graph = get_drought_module("txt/ppi.tsv",preds)

        # Mark nodes
        for node in manualseeds:
            try:
                graph.nodes[node]['isSeed'] = True
            except KeyError:
                None

        # Change node labels to uniprot names
        for node in graph.nodes:
            graph.nodes[node]['label'] = stringidconvert([node],aliasdict,'Uniprot')[0]

        nx.write_gexf(graph,os.path.join("outputs/results",file))



