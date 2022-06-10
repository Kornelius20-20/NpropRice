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
alpha = [0.1]

cutoffs = [float(i) for i in range(20,96,5)]

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
                graphfile = netprop.netprop(maingraph, seedlist, aliasfile, weight[i], alpha[j], iter[k], regen=regen)

def cutoff_graph(graphname,attr,cutoffs):

    graph = nx.read_gexf(os.path.join("outputs",graphname))

    for cut in cutoffs:
        nodes = []

        for node in graph.nodes:
            if graph.nodes[node][attr] > cut: nodes.append(node)

        newgraph = nx.subgraph(graph,nodes)

        nx.write_gexf(newgraph,os.path.join(outputdir,f"{graphname[:-5]}_{cut}.gexf"))

# generate_graphs(weight,alpha,iter)

# cutoff_graph("p-w=10-a=0.1-i=40_NP.gexf",'weight',cutoffs)

def get_top_100(graphname,attr,howmany=100):
    from get_high_scores import descendingnodes

    graph = nx.read_gexf(os.path.join("outputs", graphname))

    keys, values = zip(*nx.get_node_attributes(graph,attr).items())

    top = descendingnodes(keys,values)[:howmany]

    return top

test = get_top_100("p-w=10-a=0.1-i=40_NP.gexf",'weight')

manualseeds = "txt/manual_mined_seeds.txt"
with open(manualseeds, 'r') as file:
    manualseeds = [line.strip() for line in file.readlines()]

# add nodes to candidate list
for i in manualseeds: test.append(i)


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


    nx.write_gexf(graph,"droughtmodule.gexf")

get_drought_module("txt/ppi.tsv",test)

