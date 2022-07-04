import os
import networkx as nx


stringseeds = "txt/string_seeds.txt"

with open(stringseeds, 'r') as file:
    stringseeds = [line.strip() for line in file.readlines()]



# graph = nx.read_gexf(g)

def clusters_as_david_list(graph,attr,outputfile):
    keys,values = zip(*nx.get_node_attributes(graph,attr).items())

    clustlist = [[] for i in range(int(max(values)))]
    for i in range(len(values)):
        clustlist[int(values[i])-1].append(keys[i])

    from stringInteractions2namedInteractions import stringidconvert

    clustlist = [stringidconvert(i) for i in clustlist]

    from cluster_drought_module_greedy import transpose_lists

    clustlist = transpose_lists(clustlist)
    clustlist.insert(0,[f"cluster_{i+1}" for i in range(int(max(values)))])

    with open(outputfile,'w') as file:
        for line in clustlist:
            file.writelines("\t".join(line))
            file.writelines('\n')


def hishigaki_single_function_predictor(graph, seed_list, n=3) -> dict:
    """
        Function to return the likely candidate protein for a certain function
        (the one in the seed list) using the Hishigaki algorithm when given
        a seed protein list using the available graph
        :param graph: Built graph object of protein-protein interactions
        :param seed_list: List of seed proteins
        :param n: number of candidates to output
        :return: a dict of candidate proteins with their scores
        """


    # Find neighbors of a protein
    testprots = set(graph.nodes).difference(seed_list)

    # variables to calculate ei
    total_nodes = len(graph.nodes)
    nodes_of_function = len(set(graph.nodes).intersection(seed_list))

    candidates = {}

    for testprot in testprots:
        # Calculating ni and ei
        total_neighbors = len(list(graph.neighbors(testprot)))
        ni = len(set(graph.neighbors(testprot)).intersection(seed_list))
        ei = float(nodes_of_function * total_neighbors / total_nodes)

        try:
            x2 = ((ni - ei) ** 2) / ei  # Hishigaki equation
        # If there are no neighbors then x2 will return an error. If so move on
        # to the next protein
        except ZeroDivisionError:
            continue

        candidates[testprot] = x2

    # Code to select n best candidates only to output
    cands = {}
    for i in range(n):
        max_prot = max(candidates, key=candidates.get)
        cands[max_prot] = candidates[max_prot]
        candidates.pop(max_prot)

    return cands

def remove_strays(graph):
    # Get largest connected component
    comp = sorted(nx.connected_components(graph), key=len, reverse=True)[0]

    # Get the subgraph consisting of only the largest component
    graph = nx.subgraph(graph,comp)

    return graph

def convert_labels(graph,source = 'Uniprot'):
    from stringInteractions2namedInteractions import stringidconvert,create_aliasdict,aliasfile

    aliasdict = create_aliasdict(aliasfile)

    # Change node labels to uniprot names
    for node in graph.nodes:
        graph.nodes[node]['label'] = stringidconvert([node], aliasdict, source)[0]

    return graph


resultpath = "outputs/results"
graph = nx.read_gexf('outputs/results/p-w=100-a=0.1-i=5_NP_100.gexf')
graph = convert_labels(graph)
nx.write_gexf(graph, 'outputs/results/p-w=100-a=0.1-i=5_NP_100.gexf')

nodeseeds = {}
for node in graph.nodes:
    key = node
    try:
        value = graph.nodes[node]['isSeed']
    except KeyError:
        value = False
    nodeseeds[key] = value

from wrapper import get_top
from stringInteractions2namedInteractions import stringidconvert,aliasfile,create_aliasdict,id2stringdict

aliasdict = create_aliasdict(aliasfile)

# Get top component of graph
maincomp = remove_strays(graph)

nodes = list(get_top(maincomp, f"louvain_PC", howmany=10, returnscores=True))


formatted = []
neighbors = {}
for name, score in zip(nodes[0], nodes[1]):

    nbs = [node for node in nx.neighbors(graph,name)]
    nbmods = set([str(graph.nodes[node]['louvain']) for node in nbs])
    neighbors[name] = nbs
    formatted.append([ stringidconvert([name],aliasdict)[0] , score , nodeseeds[name], " ".join(nbmods) ])

# with open(os.path.join(resultpath, "top_hubs.csv"), 'w') as file:
#     file.write("Candidate,Score,is a seed, connecting sub-modules\n")
#
#     for name, score,seed,mods in formatted:
#         file.write(f"{name},{score},{seed},{mods}\n")

#
# path = "outputs/results"
# for _,_,files in os.walk(path):
#     for file in files:
#         if file[-5:] == ".gexf":
#             graph = nx.read_gexf(os.path.join(path,file))
#
#             clusters_as_david_list(graph, 'louvain', os.path.join(path, f"DAVID.tsv"))

newnbs = {}
for key,value in neighbors.items():
    newnbs[stringidconvert([key],aliasdict)[0]] = stringidconvert(value,aliasdict)

# Get Id converting dict
iddict = id2stringdict()

# get annotations
import gzip
annotdict = {}
with gzip.open("gz/39947.protein.info.v11.5.txt.gz",'rt') as file:
    data = [line.strip().split('\t') for line in file.readlines()]
    for line in data:
        annotdict[line[0]] = line[-1]

with open(os.path.join(resultpath,"connectedprots.tsv"),'w') as file:
    file.write("Predicted hub\tConnected nodes\tAnnotation\n")
    for key,value in newnbs.items():
        for val in value:
            file.write(f"{key}\t{val}\t{annotdict[iddict[val]]}\n")