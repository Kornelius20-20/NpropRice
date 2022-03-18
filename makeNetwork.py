"""

Author: Janith Weeraman

Code that creates the graph file from a tsv of ppi interactions and writes it into a gexf file

"""

import networkx as nx
netfile = r"txt/ppi.tsv"
seedfile = "txt/manual_mined_seeds.txt"
do_search = False
graph = nx.Graph()

def do_community_search(graph,seedfile):
    """
    Method that takes in a graph and a seedfile and does greedy community search on the graph.
    It then returns the subgraph of the original graph that corresponds to the community that contains
    the majority of the seeds in the seedfile provided

    :param graph: input graph file
    :param seedfile: as a string corresponding to the file path
    :return: graph: subgraph containing seeds
    """

    greedycom = nx.algorithms.community.greedy_modularity_communities(graph)
    # Get seed list as set
    with open(seedfile,'r') as file:
        seeds = [line.rstrip() for line in file.readlines()]
    seeds = set(seeds)

    # Find community with highest number of seeds
    results = []
    for i in range(len(greedycom)):
        results.append(len(seeds.intersection(set(greedycom[i]))))

    bestcom = results.index(max(results))

    graph = graph.subgraph(greedycom[bestcom])

    return graph

if __name__ == "__main__":

    # Open edge data file
    with open(netfile, 'r') as file:
        lines = [line.rstrip().split('\t') for line in file.readlines()]

    # create graph
    for data in lines[1:]:
        try:
            graph.add_edge(data[0], data[1], weight=data[-1])
        except IndexError:
            print(data)
            print("Something went wrong with your data file")
            break

    if do_search: graph = do_community_search(graph,seedfile)
    # write graph as gexf file
    nx.write_gexf(graph,"graph.gexf")
