"""
Attempt to code network propagation
"""
import sys

import networkx as nx
import numpy as np
import cupy as cp
import os
from stringInteractions2namedInteractions import create_aliasdict

seedlist = r"txt/string_seeds.txt"
aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
regen = False

# Load graph
graph = nx.read_gexf('graph.gexf')
# number of iterations to run to the random walk
iter = 100
# following parameters should be given in list form
# weight to give to seeds
weight = [50,100,250,500]
# restart parameter
alpha = [0.1,0.5,1.0,2.0]

# output numpy file name (without the .npy)
outfile = "p"

scale_limit = 100.0
cutoff = 50.0

def weights_from_seeds(graph,seedlist, weight=100):
    p0 = []
    # Get seed proteins and assign to them an input weight of 100 and 0 to the rest
    with open(os.path.join(seedlist), 'r') as file:
        seeds = [line.rstrip() for line in file.readlines()]

    for name in graph.nodes:
        if name in seeds:
            p0.append(100)
        else:
            p0.append(0)

    p0 = np.array(p0)

    return p0


def rwr(graph, seedlist, weight, alpha, A, invD, iter,):
    # Input weights for network
    p0 = weights_from_seeds(graph, seedlist, weight)

    invD = cp.array(invD)
    cA = cp.array(A.toarray())

    W = cp.matmul(cA, invD)

    # RWR
    alp0 = cp.array(alpha * p0)
    p = alp0
    for i in range(iter):
        p = alp0 + (1 - alpha) * cp.matmul(W, p)

    return cp.asnumpy(p)


def graph_with_weights(wgraph,alias_key,p,outcode):
    # Change graph labels to uniprot names
    p *= scale_limit / p.max()

    sub_nodes = []
    i = 0
    for node in wgraph.nodes: # for each node
        # change its label
        wgraph.nodes[node]['label'] = alias_key[node].get("Uniprot", alias_key[node]["Uniprot"])
        # add it's weight
        wgraph.nodes[node]['weight'] = p[i]

        # Get nodes to keep
        if p[i] > cutoff: sub_nodes.append(node)

        i += 1

    prunedgraph = wgraph.subgraph(sub_nodes)

    print("subgraph made")


    outgraph = outcode + '.gexf'
    nx.write_gexf(prunedgraph,outgraph)

# Get adjacency matrix of graph
A = nx.graphmatrix.adjacency_matrix(graph)

# Create degree matrix of graph
D = np.eye(graph.number_of_nodes())
i = 0
for name,degree in graph.degree:
    D[i,i] = degree
    i += 1

# Create the inverse degree matrix and save as a temp file. If a file already exists and regen = False
# then just use that instead
try:
    if regen: raise FileNotFoundError # if regen is set to True, force creation of file
    invD = np.load('tmp.npy')
except FileNotFoundError:
    invD = np.linalg.inv(D)
    with open('tmp.npy', 'wb') as file:
        np.save(file, invD)

alias_key = create_aliasdict(aliasfile)

for i in range(len(alpha)):
    # Do a random walk for some iterations and get the final weight vector for nodes
    p = rwr(graph,seedlist,weight[i],alpha[i], A,invD,iter)

    # save output numpy array
    outcode = f"{outfile} -{weight[i]}-{alpha[i]}"
    outnpy = outcode + ".npy"
    with open(outnpy,'wb') as file:
        np.save(file,p)

    graph_with_weights(graph,alias_key,p,outcode)