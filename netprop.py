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

# output numpy file name (without the .npy)
outfile = "p"

scale_limit = 100.0
cutoff = 35.0

def weights_from_seeds(graph,seedlist, weight=100):
    """
    Takes an input graph and creates a list of weights for the graph nodes, where each seed node from a list
    of seeds will get assigned the starting weight value (defaults to 100) and all other nodes will start at
    0

    :param graph: input graph
    :param seedlist: list of seed nodes to assign initial weight to
    :param weight: starting weight value
    :return: List of weights corresponding to each node in the graph
    """
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


def _rwr(p0, alpha, A, invD, iter):
    """
    Random walk with restart implementation that takes in initialization parameters and only carries out
    the walk. The initialization parameters have been taken out of the method to allow relatively constant elements
    of the calculation to be provided once without the need to recalculate multiple times. Implements cupy to
    speed up the walk steps

    :param p0: initial weights
    :param alpha: learning rate
    :param A: adjacency matrix of graph
    :param invD: inverse of the degree matrix of the graph
    :param iter: number of iterations to restart the walk
    :return: numpy array of final weights of nodes
    """

    invD = cp.array(invD)
    cA = cp.array(A.toarray())

    W = cp.matmul(cA, invD)

    # RWR
    alp0 = cp.array(alpha * p0)
    p = alp0
    for i in range(iter):
        p = alp0 + (1 - alpha) * cp.matmul(W, p)

    return cp.asnumpy(p)


def graph_with_weights(wgraph,alias_key,p,outcode='outputgraph',scale=True,cutoff=0):
    """
    Assigns node label in graph to its uniprot name, assigns its weight from the network propagation output (p)
    optionally will also scale the weights and skip labeling and adding weights to nodes that do not have a
    weight score higher than the cutoff.
    Writes the subgraph from the selected nodes into a gexf file and returns a string with the relative path
    of the written graph file

    :param wgraph: input graph
    :param alias_key: dict with keys as string IDs and values as alias names
    :param p: weights to be added to the nodes
    :param outcode: name of the output gexf file
    :param scale: whether to scale the weights. Default = True
    :param cutoff: whether to implement a cutoff. Default = 0 ie no cutoff
    :return: path to output gexf file
    """

    # Rescale the weights between 0 and the scale_limit value
    if scale:   p *= scale_limit / p.max()

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

    # Create directories to hold output files
    if not os.path.exists("/outputs"):
        os.makedirs('/outputs')
        if not os.path.exists("outputs/graphs"):
            os.makedirs('outputs/graphs')

    # write the output graph to disk
    outgraph = "outputs/graphs/" + outcode + '.gexf'
    nx.write_gexf(prunedgraph,outgraph)

    return outgraph


def netprop(graph,seedlist,aliasfile,weight,alpha,iter,scale=True,cutoff=cutoff,regen=True):


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

    createdgraphs = []
    for i in range(len(weight)):
        for j in range(len(alpha)):
            for k in range(len(iter)):
                # Input weights for network
                p0 = weights_from_seeds(graph, seedlist, weight[i])
                # Do a random walk for some iterations and get the final weight vector for nodes
                p = _rwr(p0, alpha[j], A, invD, iter[k])

                # save output numpy array
                outcode = f"{outfile}-{weight[i]}-{alpha[j]}-{iter[k]}"
                outnpy = outcode + ".npy"
                with open(outnpy,'wb') as file:
                    np.save(file,p)

                # Run the graph_with_weights and return the output filename. Append it to the list
                createdgraphs.append( graph_with_weights(graph,alias_key,p,outcode,scale,cutoff) )

    return createdgraphs
