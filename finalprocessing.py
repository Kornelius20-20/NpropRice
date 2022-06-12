import os
import networkx as nx


stringseeds = "txt/string_seeds.txt"

with open(stringseeds, 'r') as file:
    stringseeds = [line.strip() for line in file.readlines()]

manualseeds = "txt/string_seeds.txt"
with open(manualseeds, 'r') as file:
    manualseeds = [line.strip() for line in file.readlines()]


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
    clustlist.insert(0,[f"cluster_{i}" for i in range(int(max(values)))])

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
graph = nx.read_gexf('outputs/results/p-w=100-a=0.1-i=5_NP_250.gexf')
graph = convert_labels(graph)
nx.write_gexf(graph, 'outputs/results/p-w=100-a=0.1-i=5_NP_250.gexf')

nodeseeds = {}
for node in graph.nodes:
    key = graph.nodes[node]['label']
    try:
        value = graph.nodes[node]['seed']
    except KeyError:
        value = False
    nodeseeds[key] = value

from wrapper import get_top

# Get top component of graph
maincomp = remove_strays(graph)

nodes = get_top(maincomp, f"label_propagation_PC", howmany=10, returnscores=True)

nodes = list(zip(nodes[0], nodes[1]))

with open(os.path.join(resultpath, "top_hubs.csv"), 'w') as file:
    file.write("Candidate,Score,is a seed\n")
    for name, score in nodes:
        file.write(f"{name},{score},{nodeseeds[name]}\n")

print(" or ".join([i[0] for i in nodes]))


path = "outputs/results"
for _,_,files in os.walk(path):
    for file in files:
        if file[-5:] == ".gexf":
            graph = nx.read_gexf(os.path.join(path,file))

            clusters_as_david_list(graph, 'label_propagation', os.path.join(path, f"DAVID.tsv"))