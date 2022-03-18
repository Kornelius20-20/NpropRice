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



