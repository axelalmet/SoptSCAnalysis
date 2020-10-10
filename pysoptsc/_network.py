""" 
Functions that we use to construct the necessary tools for network analysis of the probability matrices
constructed using SoptSC's methods. These functions mainly rely on the NetworkX package.
"""

# Import the necessary library, i.e. NetworkX
import networkx as nx

def construct_network_from_probability_matrix(probability_matrix, weight_label = 'probability'):
    """ Constructs a directed, weighted graph from a given probability matrix.
    
    Parameters 
    ----------
    probability_matrix (sparse matrix)
        An n_cell x n_cell sparse matrix that details the individual probabilities of 
        cell-cell communication.
    edge_weight (string)
        The name of the edge weights that we would like to prescribe.

    Returns
    -------
    network (NetworkX DiGraph)
        A directed graph with n_cell nodes, which are weighted by P_ij.
    """

    # Make sure the graph is directed
    network = nx.from_scipy_sparse_matrix(probability_matrix, create_using=nx.DiGraph, edge_attribute = weight_label) 

    return network

def construct_network_list_from_probability_list(probability_matrices, weight_label = 'probability'):
    """ Constructs a directed, weighted adjacency from a given probability matrix.
    
    Parameters 
    ----------
    probability_matrices (list )
        A list of n_cell x n_cell sparse matrices that describe the individual probabilities of 
        cell-cell communication.
    edge_weight (string)
        The name of the edge weights that we would like to prescribe.

    Returns
    -------
    networks (list)
        A list of NetworkX DiGraphs with n_cell nodes, which are weighted by P_ij.
    """

    # Calculate the network for each provided probability matrix
    networks = [construct_network_from_probability_matrix(matrix, weight_label) for matrix in probability_matrices]

    return networks
