import networkx as nx
import os

outputdir = "outputs/graphs/"
algorithms = ['greedy_clusters','label_propagation','louvain']
cuts = [i for i in range(25,75,5)]


from cluster_drought_module_greedy import get_best_scoring_nodes,transpose_lists
from stringInteractions2namedInteractions import stringidconvert

def get_cands_over_cutoff(graph,attrs,cutoff):
    # For each cutoff score, for each algorithm, find which nodes score above that threshold in the most number of
    # algorithms

    algores = []
    for algo in attrs:
        res = get_best_scoring_nodes(graph,f"{algo}_PC",cutoff)
        algores.append(res)

    commons = set(algores[0])
    for ls in algores[1:]:
        commons = commons.intersection(set(ls))

    commons = list(commons)

    return commons

def add_cutoff_nodes(graph,commons,cut):
    for prot in commons:
        graph.nodes[prot][f'{cut}_candidate'] = True

for _,_,files in os.walk(os.path.join(outputdir)):

    for file in files:
        graph = nx.read_gexf(os.path.join(outputdir,file))

        commonprots = []
        for cut in cuts:
            commons = get_cands_over_cutoff(graph,algorithms,cut)
            add_cutoff_nodes(graph,commons,cut)
            commons.sort()
            commonprots.append(commons)

        with open(os.path.join("outputs/results",f"candsfor{file[:-5]}.csv"),'w') as outfile:
            lines = transpose_lists(commonprots)
            lines.insert(0,[str(i) for i in cuts])
            for line in lines:
                outfile.writelines(','.join(line))
                outfile.write('\n')
        nx.write_gexf(graph,os.path.join(outputdir,file))




def get_cand_scores(graph,candidates,attrs,output='canddata.csv'):

    outputs = []
    for algo in attrs:
        data = []
        for node in candidates:
            data.append(graph.nodes[node][f"{algo}_PC"])
        outputs.append(data)

    outputs.insert(0,candidates)

    from cluster_drought_module_greedy import transpose_lists

    outputs = transpose_lists(outputs)
    attrs.insert(0,'candidate')

    with open(output,'w') as file:
        file.writelines(','.join(attrs)+'\n')
        for line in outputs:
            line = [str(i) for i in line]
            file.writelines(','.join(line)+'\n')

get_cand_scores(graph,commons,algorithms)


def tsv_gene_clusters(graph,attr):
    clustering = nx.get_node_attributes(graph,attr)

    clusterlist = []
    for i in range(max(clustering.values())):
        clusterlist.append([key for key,value in clustering.items() if value == i+1])

    for i in range(len(clusterlist)):
        clusterlist[i] = stringidconvert(clusterlist[i])

    clusterlist = transpose_lists(clusterlist)
    with open(os.path.join(outputdir,'results',f"{attr}.txt"),'w') as file:
        # adding header with cluster numbers
        file.writelines('\t'.join([f"cluster {i+1}" for i in range(max(clustering.values()))])+'\n')

        for line in clusterlist:
            file.writelines('\t'.join(line)+'\n')

