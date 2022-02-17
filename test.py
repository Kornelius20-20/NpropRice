import networkx as nx
import numpy as np

# Get highest scoring nodes from the gexfs that aren't seeds
graph = nx.read_gexf("testgraph.gexf")


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


nodelist,weightlist = get_non_seeds(graph,True)

# Create sorted arrays of nodes and weight and print them in descending order
weightnp = np.array(weightlist)
nodelist = [nodelist[i] for i in np.argsort(weightnp)]
weightlist = [weightlist[i] for i in np.argsort(weightnp)]

print(list(reversed(nodelist)))
print(list(reversed(weightlist)))
print(len(nodelist))