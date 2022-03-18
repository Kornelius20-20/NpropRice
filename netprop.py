"""
Attempt to code network propagation
"""
import sys

import networkx as nx
import numpy as np
import cupy as cp
import os
from stringInteractions2namedInteractions import create_aliasdict

# output numpy file name (without the .npy)
outfile = "p"

scale_limit = 100.0
cutoff = 35.0

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

        # ONLY KEEPS LABELS ABOVE CUTOFF
        # Get nodes to keep
        if p[i] > cutoff: sub_nodes.append(node)

        i += 1

    prunedgraph = wgraph.subgraph(sub_nodes)

    outgraph = outcode + '.gexf'
    nx.write_gexf(prunedgraph,outgraph)

    return outgraph

def netprop(graph,seedlist,aliasfile,weight,alpha,iter,regen=False):
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
        # with open('tmp.npy', 'wb') as file:
        #     np.save(file, invD)

    alias_key = create_aliasdict(aliasfile)
    graphnames = []
    for i in range(len(weight)):
        for j in range(len(alpha)):
            for k in range(len(iter)):
                # Do a random walk for some iterations and get the final weight vector for nodes
                p = rwr(graph,seedlist,weight[i],alpha[j], A,invD,iter[k])

                # save output numpy array
                outcode = f"outputs/{outfile}-{weight[i]}-{alpha[j]}-{iter[k]}"
                outnpy = outcode + ".npy"
                with open(outnpy,'wb') as file:
                    np.save(file,p)

                graphnames.append(graph_with_weights(graph,alias_key,p,outcode))

    return graphnames