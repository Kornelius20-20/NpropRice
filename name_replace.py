import os.path, gzip

# path to ppi file
path = os.getcwd()
pinfo = "gz/39947.protein.info.v11.5.txt.gz"  # ppi file
seeds = "txt/seed_id.txt"
useSeeds = True


def getNodes(interactions, seed_nodes, depth=1):
    """
    Get the nodes that will be used to construct the drought module.
    This is a dumb method for getting all nodes of a certain depth from an initial seed list

    :param interactions:
    :param seed_nodes:
    :param depth:
    :return:
    """
    seed_nodes = set(seed_nodes)

    newnodes = []
    for i in range(depth):
        for interaction in interactions:
            # For each interaction, check which proteins are already seeds, if one of the proteins
            # is not a seed then add that protein to the new node list
            left = interaction[0] in seed_nodes
            right = interaction[1] in seed_nodes

            if left and right:
                None
            elif left:
                newnodes.append(interaction[1])
            elif right:
                newnodes.append(interaction[0])

        seed_nodes = seed_nodes.union(set(newnodes))

    return seed_nodes

with gzip.open(os.path.join(path, pinfo), 'rt') as file:
    text = [line.rstrip().split('\t') for line in file.readlines()]

protein_key = {}
for line in text:
    protein_key[line[0]] = [line[1], line[-1]]

"""
Take the new_ppi interactions and convert their string IDs into their gene names

"""
with open(os.path.join(path, "txt/new_ppi.txt"), 'r') as file:
    ppidata = [line.rstrip().split(' ') for line in file.readlines()]

if useSeeds:
    with open(os.path.join(path, seeds), 'r') as file:
        seedlist = set([line.rstrip().split()[0] for line in file.readlines()])
    seedlist = getNodes(ppidata,seedlist,1)

with open(os.path.join(path, "txt/converted_ppi.txt"), 'w') as file:
    if useSeeds:
        for line in ppidata:
            if line[0] in seedlist or line[1] in seedlist:
                file.writelines(' '.join(line) + '\n')
    else:
        for line in ppidata:
            file.writelines(' '.join(line) + '\n')