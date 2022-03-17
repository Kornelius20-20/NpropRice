import networkx as nx
import numpy as np

# Get highest scoring nodes from the gexfs that aren't seeds


def get_non_seeds(graph,withscore=True):
    # Function that

    nodelist = []
    weightlist = []

    for node,attrs in graph.nodes(data=True):
        try:
            if graph.nodes[node]['isSeed']:
                continue
        except KeyError:
            nodelist.append(graph.nodes[node]['label'])
            weightlist.append(graph.nodes[node]['weight'])
    return nodelist,weightlist

def descendingnodes(nodes,weights,returnboth=False):
    # Create sorted arrays of nodes and weight and print them in descending order
    weightnp = np.argsort(np.array(weights))
    nodelist = [nodes[i] for i in weightnp]

    if returnboth:
        weightlist = [weights[i] for i in weightnp]
        return list(reversed(nodelist)),list(reversed(weightlist))
    else: list(reversed(nodelist))


# import pandas as pd
# dframe = pd.read_csv("txt/uniprot_original.csv")
#
# for graphfile in graphfiles:
#     graphfile = graphfile[:-5] + "seedsAndGO.gexf"
#     graph = nx.read_gexf(graphfile)
#
#     # Get nodes that aren't seeds and their weights
#     nodelist,weightlist = get_non_seeds(graph,True)
#
#     # sort them in ascending order
#     nodes,weights = descendingnodes(nodelist,weightlist,True)
#
#     entries = [dframe.loc[dframe['Entry'] == node] for node in nodes]
#     results = pd.concat(entries)
#     results.to_csv(f'results-{graphfile}.csv')

