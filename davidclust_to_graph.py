import stringInteractions2namedInteractions as s2n
import networkx as nx
import os

outputdir = 'outputs/'
graphfile = 'p-w=10-a=0.5-i=75-c=35.0.gexf'
aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"

with open('geneClusterReport.txt','r') as file:
    lines = [line.rstrip().split('\t') for line in file.readlines()]

clustlines = []
# Find the lines that correspond to a start of a new cluster
for i in range(len(lines)):
    # print(lines[i])
    if lines[i][0][5:12] == "Cluster":
        clustlines.append(i)

# Generate dict of sources and their corresponding STRING IDs
sourcedict = s2n.id2stringdict(aliasfile)

graph = nx.read_gexf(os.path.join(outputdir,'graphs',graphfile))


for i in range(len(clustlines)-1):
    # Select all lines that correspond to a cluster ignoring the header line for the proteins
    clustnum = lines[clustlines[i]][0][-2:]
    clust = lines[clustlines[i]+2:clustlines[i+1]]
    for item in clust:
        item[0] = sourcedict[item[0]]
        graph.nodes[item[0]]['DAVID_clust'] = clustnum
        graph.nodes[item[0]]['DAVID_name'] = item[1]


nx.write_gexf(graph,os.path.join(outputdir,'results',graphfile))