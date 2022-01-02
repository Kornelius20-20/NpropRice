import networkx as nx
import pandas as pd
import re

graphfile = "drought_module.gexf"
dframe = "txt/processed_uniprot.csv"
interested_go = {"GO:0009415", 'GO:0009270', 'GO:0009413', 'GO:0009414', 'GO:0030912', 'GO:0051599', 'GO:0071462',
                 "GO:0009738", "GO:0009787", "GO:1901527", "GO:0009788", "GO:0009789" }

graph = nx.read_gexf(graphfile)
frame = pd.read_csv(dframe)

# Add gene ontology IDs and reviewed status as attributes
nodes = {}
for prot in graph.nodes:
        status = frame[frame['Entry'] == prot]
        if status.empty:
            None
        else:
            attributes = {'status':status['Status'].to_string()[-8:]}

            # Get GO IDs and remove the row number from the first index
            go_process = status["Gene ontology (biological process)"].to_string() \
                         + status["Gene ontology IDs"].to_string()
            # finding the GO IDs via regex search
            regex_results = re.findall("GO:\d\d\d\d\d\d\d", go_process)

            for item in regex_results:
                if item in interested_go:
                    attributes['GO'] = item
                    break

            # response to osmotic stress GO:0006970
            nodes[prot] = attributes
nx.set_node_attributes(graph, nodes)


# write graph as gexf file
nx.write_gexf(graph, "drought_module2.gexf")