"""

Author: Janith Weeraman
Date: 16/02/2021
Version: 1.6

A package with classes and methods to predict the functions
of unknown proteins in a protein-protein interaction network
using networkX

dependancies: NetworkX

"""
import networkx as nx
import gzip


def parse_seed_list(filename, delimitter='\t', header=False, protein_code_column=0) -> set:
    """
    Method to take a seed protein list as a tsv and convert only the
    protein codes into a list.

    By default the method assumes the protein codes are stored in the first
    column (index = 0)

    :param filename: Name of seed list tsv file
    :param delimitter: default is tabs
    :param header: default is False
    :param protein_code_column: the column of the tsv where protein code is
    :return: List of protein codes
    """

    with open(filename, 'r') as file:
        data = file.readlines()
        if header: data = data[1:]  # Remove header if present

    seed_list = []  # list to hold seed proteins

    for line in data:
        cols = line.strip().split(delimitter)  # split each record to columns
        seed_list.append(cols[protein_code_column])  # Add the protein to the seed list

    return set(seed_list)


def _multi_function_predictor(func: "The function to use as the predictor",
                              folder: "Folder path containing function seeds") -> dict:
    """
        Method to iteratively run the given single_function predictor on multiple seed protein files
        in a single folder. This will output a dictionary containing the likely candidate genes
        for each function (represented by each protein seed file)

        :param func: the function used as the predictor
        :param graph: a built protein-protein interaction graph
        :param folder: String path of folder where seed files are
        :return: dictionary with function (seed list) names as keys and candidates as values
        """
    import os

    # Get only the text files in the target folder
    files = os.listdir(folder)
    files = [file for file in files if file.split('.')[1] == 'txt']

    cand_functions = {}  # Dictionary to hold candidate proteins for functions

    # Iterate through file list and run predictor on each file
    for file in files:
        filename = os.path.join(folder, file)
        candidate = func(filename)  # Calling the assigned function with the graph and seed value

        cand_functions[file.split('.')[0]] = candidate

    # Returns the dictionary of functions along with their candidate genes
    return cand_functions


def make_graph(network_file, delimitter='\t'):
    """
    A method to create and return a NetworkX graph object
    created from a tsv representation of an interaction
    network.

    :param network_file: a file in tsv format of network interactions
    :param delimitter: default '\t'
    :return: a networkx graph object
    """

    with open(network_file, 'r') as file:
        lines = file.readlines()

    commented_lines = []  # variable to hold indices of commented lines

    for i in range(len(lines)):
        if lines[i][0] == '#':  # Marks commented lines
            commented_lines.append(i)
        else:
            # Splits non-commented lines into list of values
            lines[i] = lines[i].strip().split(delimitter)

    # Remove the commented lines from the list
    for i in commented_lines:
        lines.pop(i)

    # Create new graph object
    graph = nx.Graph()

    # Add edges to graph
    for row in lines:
        graph.add_edge(row[0], row[1], weight=row[-1])

    return graph


def preprocessor(networkfile, infofile) -> list:
    """
    Method that takes a download gz file for a graph and its associated protein info file to make a  set of data that
    can be used to make a graph with nodes named after the protein codes

    :param networkfile: gz file containing protein network
    :param infofile: protein info file from the string database
    :return: list where each item is a list containing the nodes and weight of a single edge in the graph
    """
    # Get the protein codes and names (without the header)
    with open(infofile, 'r') as file:
        lines = file.readlines()[1:]

    # Create a dict with codes as keys and corresponding names as values
    protein_dict = {}
    for line in lines:
        code, name = line.strip().split('\t')[:2]
        protein_dict[code] = name

    # Read the organism graph file ignoring the header
    with gzip.open(networkfile, 'rb') as file:
        protdata = file.readlines()[1:]

    # Variable for storing list of graph data
    graphdata = []

    # For each data line
    for line in protdata:
        # Convert each line from byte to text using utf-8 encoding
        line = line.decode('utf-8')
        # Get each node in the line and convert it to protein code present in infofile
        data = line.strip().split(' ')
        data[0] = protein_dict[data[0]]
        data[1] = protein_dict[data[1]]

        # Add to list of graph data
        graphdata.append(data)

    return graphdata


class NetCandPred:
    """
    Class to hold the formed graph file and containing methods to run a search for candidate proteins/genes.
    This implementation can carry out a majority voting algorithm implementation or a hishigaki algorithms
    implementation on one or many seed protein lists (if many, then the input must be the folder containing the lists)
    There is also methods to predict using each of the above algorithms, which function fits an unknown protein best

    """
    graph: nx.Graph

    def __init__(self, networkfile, infofile=None,delimitter='\t'):
        """
        Alternate constructor that can take a protein gz file and info file download straight from string-db.org to
        construct a graph directly, which can then be used with the other methods

        :param networkfile:
        :param infofile:
        """
        if infofile is None:
            self.graph = make_graph(networkfile, delimitter)
        else:
            self.graph = nx.Graph()
            graphdata = preprocessor(networkfile, infofile)
            for line in graphdata:
                self.graph.add_edge(line[0], line[1], weight=line[-1])

    def proteins_interactions(self) -> tuple:
        """
        Method to return the number of proteins and the number
        of interactions in a graph

        :param graph: graph object to query
        :return: tuple where first value is the number of proteins and
        second value is the number of interactions
        """

        return (len(self.graph.nodes), len(self.graph.edges))

    def hishigaki_single_function_predictor(self, seed_list, n=3) -> dict:
        """
            Function to return the likely candidate protein for a certain function
            (the one in the seed list) using the Hishigaki algorithm when given
            a seed protein list using the available graph

            :param graph: Built graph object of protein-protein interactions
            :param seed_list: List of seed proteins
            :param n: number of candidates to output
            :return: a dict of candidate proteins with their scores
            """

        graph = self.graph
        seed_list = parse_seed_list(seed_list)

        # Find neighbors of a protein
        testprots = set(graph.nodes).difference(seed_list)

        # variables to calculate ei
        total_nodes = len(graph.nodes)
        nodes_of_function = len(set(graph.nodes).intersection(seed_list))

        candidates = {}

        for testprot in testprots:
            # Calculating ni and ei
            total_neighbors = len(list(graph.neighbors(testprot)))
            ni = len(set(graph.neighbors(testprot)).intersection(seed_list))
            ei = float(nodes_of_function * total_neighbors / total_nodes)

            try:
                x2 = ((ni - ei) ** 2) / ei  # Hishigaki equation
            # If there are no neighbors then x2 will return an error. If so move on
            # to the next protein
            except ZeroDivisionError:
                continue

            candidates[testprot] = x2

        # Code to select n best candidates only to output
        cands = {}
        for i in range(n):
            max_prot = max(candidates, key=candidates.get)
            cands[max_prot] = candidates[max_prot]
            candidates.pop(max_prot)

        return cands

    def majority_single_function_predictor(self, seed_list, n=3) -> dict:
        """
        Function to return the likely candidate protein for a certain function
        (the one in the seed list) using the majority voting algorithm when given
        a seed protein list using the available graph

        :param graph: Built graph object of protein-protein interactions
        :param seed_list: List of seed proteins
        :param n: number of candidates to output
        :return: a dict of candidate proteins with their scores
        """
        graph = self.graph
        seed_list = parse_seed_list(seed_list)

        # Get the set of all proteins in the graph
        candidates = set(graph.nodes)

        # Find proteins that are in the graph that aren't in the seed list
        unknowns = candidates.difference(seed_list)

        scores = {}

        # For each unknown protein
        for protein in unknowns:
            # Find its neighbours and convert it into a set
            neighbours = set(graph.neighbors(protein))

            # For each neighbour present in the seed list, increment score by 1
            scores[protein] = len(neighbours.intersection(seed_list))

        # Return the protein that has the highest neighbour score according to majority voting
        cands = {}
        for i in range(n):
            max_prot = max(scores, key=scores.get)
            cands[max_prot] = scores[max_prot]
            scores.pop(max_prot)

        return cands

    def majority_multi_function(self, folder: "Folder path containing function seeds") -> dict:
        """
        Method to run the majority voting algorithm on multiple seed lists in a folder and return the output

        :param folder: folder containing seed lists
        :return: dictionary of function: candidate pairs
        """
        return _multi_function_predictor(self.majority_single_function_predictor, folder)

    def hishigaki_multi_function(self, folder: "Folder path containing function seeds") -> dict:
        """
            Method to run the hishigaki algorithm on multiple seed lists in a folder and return the output

            :param folder: folder containing seed lists
            :return: dictionary of function: candidate pairs
        """
        return _multi_function_predictor(self.hishigaki_single_function_predictor, folder)

    def majority_best_function_predictor(self, protein: "unknown protein", folder: "folder of function seeds") -> list:
        """
        Method to predict the most accurate function for an unknown protein by considering all seed lists
        given in a folder one at a time using the majority voting algorithm

        :param protein: Unknown protein which to check function of
        :param folder: folder containing seed lists
        :return: list containing most likely function(s) of unknown protein
        """

        import os
        # Get only the text files in the target folder
        files = os.listdir(folder)
        files = [file for file in files if file.split('.')[1] == 'txt']

        graph = self.graph
        # Get neighbours of unknown protein
        neighbours = set(graph.neighbors(protein))

        # Find the function that has the most neighboring proteins for the unknown protein
        likely_function: str
        n_bors = 0
        for file in files:
            filename = os.path.join(folder, file)
            seeds = parse_seed_list(filename)
            num = len(neighbours.intersection(seeds))  # Get number of seed proteins that are neighbors

            if num > n_bors:
                likely_function = file.split('.')[0]  # Assigning max neighbor as likely function
                n_bors = num
            elif num == n_bors:
                likely_function = likely_function + ', ' + file.split('.')[0]  # for equal max values

        # return a list with length equal to number of functions that had max number of neighbors
        if ',' in likely_function:
            return likely_function.split(', ')
        else:
            return [likely_function]

    def hishigaki_best_function_predictor(self, protein: "unknown protein", folder: "folder of function seeds") -> list:
        """
        Method to predict the most accurate function for an unknown protein by considering all seed lists
        given in a folder one at a time using the majority voting algorithm

        :param protein: Unknown protein which to check function of
        :param folder: folder containing seed lists
        :return: list containing most likely function(s) of unknown protein
        """
        import os
        # Get only the text files in the target folder
        files = os.listdir(folder)
        files = [file for file in files if file.split('.')[1] == 'txt']

        graph = self.graph
        neighbors = set(graph.neighbors(protein))  # Get the neighbours around the protein of interest

        # variables to calculate ei
        total_nodes = len(graph.nodes)
        total_neighbors = len(neighbors)

        best_function = ''
        function_x2 = 0

        for file in files:
            filename = os.path.join(folder, file)
            seeds = parse_seed_list(filename)

            nodes_of_function = len(set(graph.nodes).intersection(seeds))
            ni = len(neighbors.intersection(seeds))

            ei = nodes_of_function * total_neighbors / total_nodes

            x2 = ((ni - ei) ** 2) / ei  # Hishigaki equation

            if x2 > function_x2:
                best_function = file.split('.')[0]
                function_x2 = x2
            elif x2 == function_x2:
                best_function = best_function + ', ' + file.split('.')[0]

        # return a list with length equal to number of functions that had max number of neighbors
        if ',' in best_function:
            return best_function.split(', ')
        else:
            return [best_function]
