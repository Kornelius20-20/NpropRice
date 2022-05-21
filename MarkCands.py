import networkx as nx
import os

outputdir = "outputs/graphs/"
algorithms = ['greedy_clusters','label_propagation','louvain']
cuts = [60]


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

        for cut in cuts:
            commons = get_cands_over_cutoff(graph,algorithms,cut)
            add_cutoff_nodes(graph,commons,cut)

        nx.write_gexf(graph,os.path.join(outputdir,file))

from test2 import get_cand_scores

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

def getGOofClusters():
    from suds.client import Client

    # Getting functional annotation chart from DAVID for a gene list and figuring out how to parse it
    genetype = "UNIPROT_ID"
    tool = "term2term"
    annotlist = ['GOTERM_BP_FAT','GOTERM_CC_FAT','GOTERM_MF_FAT','INTERPRO','KEGG_PAHWAY']
    idlist = ['Q6K9M5_ORYSJ', 'Q93X08_ORYSJ', 'Q0J0N3_ORYSJ', 'SUS3_ORYSJ', 'C7J2A4_ORYSJ', 'SUS4_ORYSJ', 'SUS5_ORYSJ',
              'SPSA3_ORYSJ', 'SPSA4_ORYSJ', 'Q6ZGL5_ORYSJ', 'SPSA5_ORYSJ', 'Q6ERD9_ORYSJ', 'Q8GTK0_ORYSJ', 'Q2MJA8_ORYSJ',
              'SPSA1_ORYSJ', 'A0A0P0V0U6_ORYSJ', 'Q75II7_ORYSJ', 'Q6ZCH3_ORYSJ', 'Q2MJA7_ORYSJ', 'Q10P27_ORYSJ',
              'Q7F5L2_ORYSJ', 'SUS7_ORYSJ', 'Q5JNJ1_ORYSJ', 'A0A0N7KG85_ORYSJ', 'Q9AX07_ORYSJ', 'Q0DTU6_ORYSJ',
              'A0A0P0UZR1_ORYSJ', 'Q6Z548_ORYSJ', 'SPSA2_ORYSJ', 'SSG1_ORYSJ', 'SUS2_ORYSJ', 'SUS1_ORYSJ', 'SUS6_ORYSJ',
              'Q0DDZ4_ORYSJ']

    annot = ','.join(annotlist)
    ids = ','.join(idlist)
    listType = 0

    serviceurl = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'

    # Create service client using WSDL
    client = Client(serviceurl)
    client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

    # authenticate
    client.service.authenticate('2017s16486@stu.cmb.ac.lk')

    client.service.addList(ids, genetype, "testlist", listType)
    result = client.service.getTermClusterReport()

    print(result)