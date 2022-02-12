"""

original go = 9414 response to water deprivation

"""
import json,re
import networkx as nx

filename = "droughtGOtree.json"

graph = nx.DiGraph()

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

            nodeid = result[0]
            if int(nodeid[-4:]) != 7610:
                # Create edges containing GO terms associated with current term
                graph.add_edge(id,nodeid)

                result = re.search("\(.*?\)", relstring)
                graph.nodes[nodeid]['label'] = result[0][1:-1]

    graph.nodes[id]['label'] = name

nx.write_gexf(graph,'GOgraph.gexf')