"""

original go = 9414 response to water deprivation

"""
import json,re
import networkx as nx

filename = "droughtGOtree.json"

graph = nx.Graph()

# load json file for a GO term as a dict
with open(filename,'r') as file:
    dictjson = json.load(file)
    test = dictjson['results']

# Create a relations graph of the GO terms
for item in test:
    id = item['id']
    name = item['name']
    definition = item['definition']

    for relation in item['history']:
        if relation['category'] == 'RELATION':
            relstring = relation['text']
            result = re.search("GO:\d\d\d\d\d\d\d",relstring) # Find GO term in relation

            # Create edges containing GO terms associated with current term
            graph.add_edge(id,result[0])

    graph.nodes[id]['label'] = name

nx.write_gexf(graph,'GOgraph.gexf')