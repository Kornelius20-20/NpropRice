import networkx as nx
from scipy.sparse import csc_array
import os


def weights_from_seeds(graph, seeds):
    """
    Takes an input graph and creates a list of weights for the graph nodes, where each seed node from a list
    of seeds will get assigned the starting weight value (defaults to 100) and all other nodes will start at
    0

    :param graph: input graph
    :param seedlist: list of seed nodes to assign initial weight to
    :param weight: starting weight value
    :return: List of weights corresponding to each node in the graph
    """

    p0 = []

    for name in graph.nodes:
        if name in seeds:
            p0.append(nx.degree(graph,name))
        else:
            p0.append(0)

    p0 = csc_array(p0)

    return p0.transpose()
def netprop(graph,seeds,alpha = 0.7,num_iters=0,threshold = 1.0e-10):
    """
    Given a networkx graph and a list of seeds, assigns a weight to each seed and propagates it out into
    the network for a certain number of iterations or until a threshold for the sum of node weights is 
    reached.

    :param graph: input graph
    :param seedlist: list of seed nodes
    :param alpha: learning rate
    :param iter: number of iterations to restart the walk
    :param threshold: limit of change of weight at which the rwr is stopped
    :return: numpy array of final weights of nodes
    
    """


    A = nx.adjacency_matrix(graph).tocsc()
    A = csc_array(A,dtype="float64")

    # Normalize adjacency matrix
    for i in range(A.shape[1]):
        colsum = A.data[A.indptr[i]:A.indptr[i + 1]].sum()

        if colsum > 0:
            A.data[A.indptr[i]:A.indptr[i + 1]] = (A.data[A.indptr[i]:A.indptr[i + 1]]/colsum)


    p0 = weights_from_seeds(graph,seeds)

    p = alpha * p0
    const_fact = p.copy()

    if num_iters > 0:
        for i in range(num_iters):
            _p = p.copy()
            p = const_fact + (1-alpha)*A.dot(p)
            print(f"change in probability: {p.sum() - _p.sum()}")

    else:
        num_iters = 1

        p = const_fact + (1-alpha)*A.dot(p)
        _p = csc_array(p.shape, dtype=float)

        while (p.sum() - _p.sum()) > threshold:
            _p = p.copy()
            p = const_fact + (1-alpha)*A.dot(p)
            num_iters += 1

            print(f"change in probability: {p.sum() - _p.sum()}")


    print(f"RWR ran for {num_iters} iterations")


    return p.toarray()
