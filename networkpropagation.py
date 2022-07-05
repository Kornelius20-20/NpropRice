"""

Author: Janith Weerman
Date: 08/06/2022

Wrapper script for carrying out network propagation

"""

import os,argparse,gzip
import networkx as nx
from netprop import netprop


def generate_graphs(maingraph, seedlist, weight, alpha, iter, regen=False,mark=True):
    for i in range(len(weight)):
        for j in range(len(alpha)):
            for k in range(len(iter)):
                # Run network propagation with the given values
                graphfile = netprop(maingraph, seedlist, weight[i], alpha[j], iter[k], regen=regen)

                if mark: graphfile = add_seeds(graphfile,seedlist)

                nx.write_gexf(graphfile,os.path.join("outputs/graphs",
                                                     f"propagated w-{weight[i]},a-{alpha[j]},i-{iter[k]}.gexf"))

def add_seeds(graph,seedlist):
    # Mark nodes
    for node in seedlist:
        try:
            graph.nodes[node]['isSeed'] = True
        except KeyError:
            None

    return graph

if __name__ == "__main__":

    # Create directories to hold output files
    if not os.path.exists("/outputs"):
        os.makedirs('/outputs')
    if not os.path.exists("outputs/graphs"):
        os.makedirs('outputs/graphs')
    if not os.path.exists("outputs/results"):
        os.makedirs('outputs/results')

    # handling command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('graph') # path to gz file containing graph data
    parser.add_argument('seeds') # a plain text list of seed proteins as STRING IDs
    parser.add_argument("-w", default=100) # weights
    parser.add_argument("-a", default=0.1) # learning parameters
    parser.add_argument("-i", default=5) # iterations
    parser.add_argument("-mark", default=True) # Mark the seeds in the final network
    parser.add_argument("-reg", default=False)  # regen
    parser.add_argument("-cut", default=0.7) # cutoff for graph generation from the gz
    args = parser.parse_args()

    # Convert each argument to list form to preserve compatibility with multiple inputs for each parameter
    args.w = [float(i) for i in args.w.split(",")] if type(args.w) == str else [args.w]
    args.a = [float(i) for i in args.a.split(",")] if type(args.a) == str else [args.a]
    args.i = [float(i) for i in args.i.split(",")] if type(args.i) == str else [args.i]
    args.reg = bool(args.reg)

    # Create graph from file
    maingraph = nx.Graph()
    # Open gzip file line by line and add it to graph
    with gzip.open(os.path.join(args.graph),'rt') as file:
        # Loop through all lines in a file and get only protein names and combined score
        newcutoff = args.cut * 1000
        while True:
            line = file.readline()
            if line:
                line = line.rstrip().split(' ')
                formatted = [line[i] for i in [0, 1, -1]]
                try:
                    # Only write files that are above the cutoff
                    if int(formatted[-1]) >= newcutoff: maingraph.add_edge(formatted[0], formatted[1], formatted[-1])
                except ValueError:
                    # Exception for the first line that contains header data
                    None

            else:
                break  # break out of the loop after last line is written



    # seed list
    with open(args.seeds,'r') as file:
        seedslist = [seed.strip() for seed in file.readlines()]

    generate_graphs( maingraph,
                     seedslist,
                     args.w,
                     args.a,
                     args.i,
                     args.reg,
                     args.mark )