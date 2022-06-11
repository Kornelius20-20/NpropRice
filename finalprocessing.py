import os
import networkx as nx


stringseeds = "txt/string_seeds.txt"
g = "outputs/p-w=10-a=0.1-i=40_NP.gexf"


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

def do_hishi():
    manualseeds = "txt/string_seeds.txt"
    with open(manualseeds, 'r') as file:
        manualseeds = [line.strip() for line in file.readlines()]

    graph = nx.read_gexf("outputs/p-w=100-a=0.1-i=5_NP.gexf")

    keys,values = zip(*hishigaki_single_function_predictor(graph,manualseeds,500).items())

    from get_high_scores import descendingnodes

    results,nums = descendingnodes(keys,values,True)

    results.extend(manualseeds)

    from wrapper import get_drought_module,add_seeds

    graph = get_drought_module("txt/ppi.tsv",results)
    add_seeds(graph,manualseeds)

    nx.write_gexf(graph,"test.gexf")

for _,_,files in os.walk("outputs/results"):

    for file in files:
        if file[-5:] == ".gexf":

            graph = nx.read_gexf(os.path.join("outputs/results",file))

            clusters_as_david_list(graph,'label_propagation_PC',os.path.join("outputs/results",f"{file[:-5]}.csv"))

def add_as_seeds_code():
    for _,_,files in os.walk("outputs"):

        for file in files:

            graph = nx.read_gexf(os.path.join("outputs",file))

            for node in manualseeds:
                try:
                    graph.nodes[node]['isSeed'] = True
                except KeyError:
                    None

            nx.write_gexf(graph,os.path.join("outputs",file))

        break