import NetCandPred,os
from networkx.algorithms.link_analysis import pagerank
network = NetCandPred.NetCandPred(os.path.join("txt/output_file.txt"),delimitter=' ')

print(pagerank(network.graph))