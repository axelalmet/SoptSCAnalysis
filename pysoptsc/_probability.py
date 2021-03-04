""" 
Functions that we use to calculate the cell-cell communication probabilities. Based off the methods
from Wang et al. (2020).
"""

import numpy as np
import pandas as pd
import itertools as it
from scipy.sparse import csr_matrix

def get_expression_matrix(ann_data, gene):
    """ Method to obtain the resulting gene expressions for specified genes

    Parameters
    ----------
    ann_data (annotated data frame)
        The annotated data matrix, which contains many of the attributes we will use to calculate the
        communication probabilities.
    genes (string)
        The gene for which we want the (normalised to 1) expressions.

    Returns
    -------
    gene_expressions (sparse matrix)
        An n_cell x n_gene expression matrix.
    """

    num_cells = ann_data.n_obs # Get the number of cells
    gene_index = ann_data.var_names.get_loc(gene) # Get the location of the gene in the matrix (column index)

    gene_expressions = csr_matrix(ann_data.X[:, gene_index]).transpose() # Expression corresponds to the relevant column of the matrix

    if np.max(gene_expressions) != 0.0:
        gene_expressions /= gene_expressions.max()

    return gene_expressions.reshape((num_cells, 1))

def get_mean_expression(ann_data, genes, mean = 'Arithmetic'):
    """ Method to obtain the resulting gene expressions for specified genes

    Parameters
    ----------
    ann_data (annotated data frame)
        The annotated data matrix, which contains many of the attributes we will use to calculate the
        communication probabilities.
    genes (list of strings)
        The list of genes for which we want the (normalised to 1) expressions.
    mean (string)
        The type of mean we will calculate. Should be either the Arithmetic or Geometric mean.

    Returns
    -------
    mean_expressions (sparse matrix)
        An n_cell x n_gene expression matrix.
    """

    # Make sure that the specified mean is correct
    if mean not in ['Arithmetic', 'Geometric']:
        raise ValueError('Specified mean type must be either Arithmetic or Geometric.')

    num_cells = ann_data.n_obs # Get the number of cells

    if type(genes) != list: # If we have only specified one gene, turn this into a list to make things easier.
        genes = [genes]

    num_genes = len(genes) # Get the number of genes

    # Initialise the mean expression vector
    mean_expressions = np.zeros((num_cells, 1))

    if mean == 'Arithmetic':

        mean_expressions = csr_matrix((num_cells, 1)) # Initialise the vector of zeros

        for gene in genes:

            gene_expressions = get_expression_matrix(ann_data, gene) # Get the corresponding gene expression vector

            mean_expressions += gene_expressions # Update the mean expressions

        return mean_expressions / num_genes # Divide by the number of genes

    else: # Should be geometric

        mean_expressions = csr_matrix(np.ones((num_cells, 1))) # Initialise the vector (should be all ones for products)

        for gene in genes:

            gene_expressions = get_expression_matrix(ann_data, gene) # Get the corresponding gene expression vector

            mean_expressions *= gene_expressions # Update the mean expressions

        return mean_expressions**(1.0 / num_genes) # Take the n_genes root

def calculate_communication_probabilities(ann_data, ligand_expressions, receptor_expressions,\
                                            upregulated_gene_expressions = None, downregulated_gene_expressions = None):
    """ 
    Method to calculate the cell-cell communication probabilities for a given ligand-receptor pair.

    Parameters
    ---------
    ann_data (annotated data frame)
        The annotated data matrix, which contains many of the attributes we will use to calculate the
        communication probabilities.
    ligand_expressions (sparse n_cell x 1 matrix)
        The normalised vector of ligand expressions that we will consider
    receptor_expressions (sparse n_cell x 1 matrix)
        The normalised vector of receptor expressions that we will consider.
    upregulated_gene_expressions (sparse n_cell x 1 matrix)
        The list of downstream genes we believe to be upregulated by the ligand-receptor pair,
        if there are any specified.
    downregulated_gene_expressions (sparse n_cell x 1 matrix)
        The list of downstream genes we believe to be downregulated by the ligand-receptor pair,
        if there are any specified.
    Returns
    ------
    communication_probabilities (sparse matrix)
        An n_cell x n_cell communication matrix P, where P_{ij} is the probability that cell i communicates
        with cell j.
    """

    # The number of cells corresponds to the number of rows
    num_cells = ann_data.n_obs # First dimension corresponds to rows
    
    # Initialise the communication matrix
    communication_probabilities = np.zeros((num_cells, num_cells))

    # Calculate the matrix for ligand-receptor interactions
    ligand_diag = np.diag(ligand_expressions.flatten()) # Diagonal matrix with the ligand expressions
    receptor_matrix = np.tile(receptor_expressions.transpose(), (num_cells, 1)) # Matrix with receptor values as the rows

    communication_probabilities = np.exp( - 1.0 * (np.dot(ligand_diag, receptor_matrix) + 1e-12)**(-1.0)) # Fudge factor of 1e-12 so we don't divide by zero

    # If the upregulated gene matrix has been specified, calculate the means in each column
    if upregulated_gene_expressions is not None:

        # Calculate the quantificaiton of the upregulated downstream genes
        beta = np.exp(-1.0 * (upregulated_gene_expressions + 1e-12)**(-1.0) )

        # Calculate the scaling of alpha due to beta
        beta_matrix = np.tile(beta.transpose(), (num_cells, 1))
        Kappa = np.dot(communication_probabilities, (communication_probabilities + beta_matrix + 1e-12)**(-1.0) ) # Fudge factors to avoid division by zero

        upstream_scaling = np.multiply(Kappa, beta_matrix) # This ensures that the values of beta multiply with the columns of Kappa

        communication_probabilities = np.multiply(communication_probabilities, upstream_scaling) # Multiply the probabilities element-wise by the scaling

    # Construct the target gene matrix of downregulated genes, if specified.
    if downregulated_gene_expressions is not None:

        # Calculate the quantification of hte upregulated downstream genes
        gamma = np.exp(-1.0 * downregulated_gene_expressions)

        # Calculate the scaling of alpha due to gamma
        gamma_matrix = np.tile(gamma.transpose(), (num_cells, 1))
        Lambda = np.dot(communication_probabilities, (communication_probabilities + gamma_matrix + 1e-12)**(-1.0) ) # Fudge factors to avoid division by zero

        downstream_scaling = np.multiply(Lambda, gamma_matrix) # This makes sure the values of gamma multiplies with the columns of Lambda

        communication_probabilities = np.multiply(communication_probabilities, downstream_scaling) # Multiply the probabilities element-wise by the scaling

    # Truncate probabilities that are too small
    communication_probabilities[communication_probabilities < 1e-6] = 0.0

    # We now divide each non-zero row by the sums
    nonzero_rows = communication_probabilities.nonzero()[0]

    communication_probabilities[nonzero_rows, :] /= communication_probabilities[nonzero_rows, :].sum(axis=1, keepdims=True) 

    return csr_matrix(communication_probabilities) # Return the matrix

def calculate_communication_probabilities_list(ann_data, ligand_receptor_pairs, upregulated_genes = [], downregulated_genes = [], mean = 'Arithmetic'):
    """ 
    Method to calculate the individual probabilities for each ligand-receptor pair
    within a signalling pathway.

    Parameters
    ---------
    ann_data (annotated data frame)
        The annotated data matrix, which contains many of the attributes we will use to calculate the
        communication probabilities.
    ligand_receptor_pairs (list of string tuples)
        The list of ligand-receptor pairs.
    upregulated_genes (list of strings)
        The list of downstream genes we believe to be upregulated by the ligand-receptor pair,
        if there are any specified.
    downregulated_genes (list of strings)
        The list of downstream genes we believe to be downregulated by the ligand-receptor pair,
        if there are any specified.
    mean ('Arithmetic' or 'Geometric)
        Specifies how we calculate the mean ligand or mean receptor expression.

    Returns
    ------
    communication_probabilities (list)
        A list of sparse n_cell x n_cell communication matrices P, where P_{ij} is the probability that cell i communicates
        with cell j for a ligand-receptor pair.
    """

    # Construct the mean ligand expressions
    ligands = [pair[0] for pair in ligand_receptor_pairs]
    receptors = [pair[1] for pair in ligand_receptor_pairs]

    mean_ligand_expressions = [get_mean_expression(ann_data, ligand, mean) for ligand in ligands]
    
    # Construct the mean receptor expressions
    mean_receptor_expressions = [get_mean_expression(ann_data, receptor, mean) for receptor in receptors]

    # Initialise the downstream gene matrices, so that we can pass them into the function later
    upregulated_gene_expressions = None
    downregulated_gene_expressions = None

    # If downstream genes have been specified, we should get these matrices too
    if upregulated_genes:
        upregulated_gene_expressions = get_mean_expression(ann_data, upregulated_genes).toarray() # Default should be arithmetic mean

    if downregulated_genes:
        downregulated_gene_expressions = get_mean_expression(ann_data, downregulated_genes).toarray() # Default should be arithmetic mean

    individual_probabilities = []

    # Calculate the communication probabilities for each ligand-receptor pair
    num_pairs = len(ligand_receptor_pairs)

    for index in range(num_pairs):

        # Calculate the corresponding probability matrix
        ligands = mean_ligand_expressions[index].toarray()
        receptors = mean_receptor_expressions[index].toarray() 
        individual_probabilities.append(calculate_communication_probabilities(ann_data,\
                                                                            ligands, receptors,\
                                                                            upregulated_gene_expressions,\
                                                                            downregulated_gene_expressions))
        
    return individual_probabilities

def calculate_aggregated_probability_matrix(individual_probability_matrices):
    """ 
    Method to calculate the aggregated probability matrix for a signalling pathway.

    Parameters
    ---------
    individual_probability_matrices (list)
        A list of the n_cell x n_cell communication matrices P_{ij}, which describe the probability
        that cell i communicates with cell j for a given ligand-receptor pair. The number of matrices
        corresponds to the number of ligand-receptor pairs within a signalling pathway.

    Returns
    ------
    aggregated_probabilities (sparse matrix)
        An n_cell x n_cell communication matrix Pagg, where Pagg_{ij} is the aggregated
        probability that cell i communicates with cell j.
    """

    num_ligand_receptor_pairs = len(individual_probability_matrices) # Length of the list should be the number of ligand-receptor pairs
    matrix_shape = individual_probability_matrices[0].shape # We're assuming they all have the same shape (they should!)

    # Initialise the aggregated probability matrix
    aggregated_probabilities = csr_matrix(matrix_shape)

    # Average over all matrices
    for matrix in individual_probability_matrices:
        aggregated_probabilities += matrix
    
    return aggregated_probabilities / num_ligand_receptor_pairs 

def normalise_aggregated_probabilities(aggregated_probability_matrix):
    """
    Method to normalise the aggregated probability matrices so that rows sum to one.
    We leave this as an optional function, just to distinguish between what is calculated in
    the original Wang et al. (2019) paper and what we have chosen to do.
    
    Parameters
    ---------
    aggregated_probability_matrix (sparse matrix)
        The aggregated probability matrix in question.

    Returns
    -------
    normalised_matrix (sparse matrix)
        The normalised equivalent of the aggregated matrix.
    """

    normalised_aggregated_probability_matrix = aggregated_probability_matrix

    # Get the non-zero row and column indices of the matrix
    nonzero_indices = aggregated_probability_matrix.nonzero()

    nonzero_rows = nonzero_indices[0] # Non-zero row indices

    normalised_aggregated_probability_matrix[nonzero_rows, :] /= normalised_aggregated_probability_matrix[nonzero_rows, :].sum(axis=1, keepdims=True) 

    return normalised_aggregated_probability_matrix

def calculate_cluster_probability_matrix(ann_data, probability_matrix, cluster_label):
    """
    Method to calculate the cluster-cluster probability matrix for a given ligand-receptor pair.

    Parameters
    ----------
    ann_data (annotated dataframe)
        The annotated dataframe that stores information like the cell names, gene IDs, cluster labels, etc.
    individual_probability_matrix (sparse matrix)
        The n_cell x n_cell communication probability matrix P_{ij}, describing the probability that
        cell i communicates with cell j. This could apply to either an individual ligand-receptor pair
        or to the aggregated probability matrix.
    cluster_label (string)
        The name of the cluster label sthat we want to consider.
    Returns
    --------
    cluster_probability_matrix (dense matrix)
        The n_cluster x n_cluster communication probability matrix P_{AB}, describing the probability that
        cluster A communicates with cluster B via the specified ligand-receptor pair.
    """
    cell_clusters = ann_data.obs[cluster_label] # Get the clusters for each cell label
    cluster_labels = list(cell_clusters.unique()) # Get the unique list of clusters
    num_clusters = len(cluster_labels)

    cluster_pairs = {(i, j):1 for i, j in it.product(range(num_clusters), range(num_clusters))} # It's faster to construct the pairs like this

    cluster_probability_matrix = np.zeros((num_clusters, num_clusters))

    # The general assumption is that there shouldn't be too many clusters, so I think we can just do a nested for loop
    # (we can speed this up if we need to via itertools and dictionaries.)
    for pair in cluster_pairs:

        # Get the clusters
        A, B = pair
        cluster_A = cluster_labels[A]
        cluster_B = cluster_labels[B]

        # Get the indices corresponding to each cluster
        cells_A = np.where(cell_clusters == cluster_A)[0]
        cells_B = np.where(cell_clusters == cluster_B)[0]

        # Get the possible porbabilities
        probabilities_A_B = np.array([probability_matrix[i, j] for i, j in it.product(cells_A, cells_B)]) # It's faster to construct the pairs like this

        cluster_probability_matrix[A, B] += probabilities_A_B.sum()

    # We now divide each non-zero row by the sums (to avoid 0/0)
    nonzero_rows = cluster_probability_matrix.nonzero()[0]

    if len(nonzero_rows) > 0:
        cluster_probability_matrix[nonzero_rows, :] /= cluster_probability_matrix[nonzero_rows, :].sum(axis=1, keepdims=True) 

    return cluster_probability_matrix