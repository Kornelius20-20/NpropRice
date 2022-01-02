import networkx as nx
import pandas as pd

graphfile = "graph.gexf"
dframe = "txt/processed_uniprot.csv"

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

greedycom = nx.community.greedy_modularity_communities(graph)

# add the community each node belongs to as a node attribute
coms = {}
for i in range(len(greedycom)):
    for prot in greedycom[i]:
        status = frame[frame['Entry'] == prot]
        if status.empty:
            graph.add_node(prot, community=i, )#status='N/A',GO='N/A')
        else:
            # Get GO IDs and remove the row number from the first index
            hold = status['Gene ontology IDs'].to_string().split(';')
            hold[0] = hold[0][-10:]
            hold = ';'.join(hold)
            graph.add_node(prot, community=i, )#status=status['Status'].to_string()[-8:],GO=hold)



# write graph as gexf file
nx.write_gexf(graph, "drought_module.gexf")
