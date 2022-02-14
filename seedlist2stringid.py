"""

Author: Janith Weeraman

method that converts a named protein list to a string Id list based on an aliasfile

"""

import os,gzip

aliases = "gz/39947.protein.aliases.v11.5.txt.gz"
seedfile = "txt/manual_mined_seeds.txt"
seedout = "txt/string_seeds.txt"


def prot2stringid(protlist, aliasfile):
    """
    method that takes a list of proteins as a seed file and returns a list of STRINGIDs corresponding to the proteins
    after looking them up against a given aliasfile

    :param protlist: list of proteins as strings
    :param aliasfile: protein alias gzip file from STRING website
    :return: list of STRINGIDs
    """
    alias_key = {}
    with gzip.open(os.path.join(aliasfile), 'rt') as file:
        line = file.readline()

        # Create a dict with Alias: Corresponding ID
        while line:
            items = line.rstrip().split()[:2]
            alias_key[items[1]] = items[0]
            line = file.readline()

    # Convert seed list to STRING ID list
    seedlist = []
    for prot in protlist:
        if prot in alias_key.keys(): seedlist.append(alias_key[prot])

    return seedlist


if __name__ == "__main__":
    with open(os.path.join(seedfile), 'r') as file:
        protlist = [prot.rstrip() for prot in file.readlines()]

    seedlist = prot2stringid(protlist,aliases)

    with open(os.path.join(seedout),'w') as file:
        file.writelines('\n'.join(seedlist))