import networkx as nx
netfile = r"txt/ppi.tsv"
seedfile = "txt/manual_mined_seeds.txt"
do_search = False
graph = nx.Graph()

# Open edge data file
with open(netfile,'r') as file:
    lines = [line.rstrip().split('\t') for line in file.readlines()]

# create graph
for data in lines[1:]:
    try:
        graph.add_edge(data[0], data[1], weight=data[-1])
    except IndexError:
        print(data)
        break


def do_community_search(graph):

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

if do_search: graph = do_community_search(graph)
# write graph as gexf file
nx.write_gexf(graph,"graph.gexf")


# THIS PART IS NO LONGER USED
# rescale the pagerank output to between 0-255

def rescale(valuedict,scalemax=256,mincutoff=0):
    """
    Method to rescale a dictionary with numeric values

    :param valuedict: the dictionary with keys and numeric values
    :param scalemax: the max value to scale to
    :param mincutoff: the value below which inputs will be rounded to 0
    :return:
    """

    scalemax = scalemax -1
    smol = min(valuedict.values())
    big = max(valuedict.values())
    for key in valuedict.keys():
        newval = int(valuedict[key]/(big-smol) * scalemax)
        valuedict[key] = newval if newval >= mincutoff else 0
