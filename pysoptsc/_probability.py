""" 
Functions that we use to calculate the cell-cell communication probabilities. Based off the methods
from Wang et al. (2020).
"""

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

def calculate_communication_probabilities(ann_data, ligand_receptor_pair, upregulated_genes = [], downregulated_genes = []):
    """ 
    Method to calculate the cell-cell communication probabilities for a given ligand-receptor pair.

    Parameters
    ---------
    ann_data (annotated data frame)
        The annotated data matrix, which contains many of the attributes we will use to calculate the
        communication probabilities.
    ligand_receptor_pair (string tuple)
        The ligand-receptor pair, in that order.
    upregulated_genes (list of strings)
        The list of downstream genes we believe to be upregulated by the ligand-receptor pair,
        if there are any specified.
    downregulated_genes (list of strings)
        The list of downstream genes we believe to be downregulated by the ligand-receptor pair,
        if there are any specified.
    Returns
    ------
    communication_probabilities (sparse matrix)
        An n_cell x n_cell communication matrix P, where P_{ij} is the probability that cell i communicates
        with cell j.
    """

    # Get the raw gene expression data
    gene_expression_data = ann_data.X 

    # The number of cells corresponds to the number of rows
    num_cells = ann_data.n_obs # First dimension corresponds to rows
    
    # Initialise the communication matrix
    communication_probabilities = csr_matrix((num_cells, num_cells))

    # We need to find which column corresponds to the ligand and receptor, respectively.
    # NB One caveat of this method currently is we're assuming that the names of the ligands, receptors
    # genes all are in the same format as what they are stored in the annotated data. That is, they're
    # either all in capitals, or only the first character is capitalised. We should probably introduce
    # something to catch this in the long term.

    gene_ids = ann_data.var_names # Get the list of gene IDs

    ligand = ligand_receptor_pair[0] # Get the name of the ligand
    ligand_column_index = gene_ids.get_loc(ligand) # Find where this ligand occurs along the columns of the matrix 
    ligand_expression_data = gene_expression_data[:, ligand_column_index] # Returns n_cell x 1 sparse matrix of ligand expression levels

    receptor = ligand_receptor_pair[1] # Get the name of the receptor
    receptor_column_index = gene_ids.get_loc(receptor) # Find the column index that corresponds to this receptor
    receptor_expression_data = gene_expression_data[:, receptor_column_index] # Returns n_cell x 1 sparse matrix of receptor expression levels\

    # We only proceed with the probability calculations if BOTH the ligand and receptor gene expressions are non-zero for a positive number of cells.
    # This is based on the underlying assumptions of the method from Wang et al. (2019)
    nonzero_ligand_indices = np.nonzero(ligand_expression_data)[0] # Get the array of non-zero indices for ligand expressions
    nonzero_receptor_indices = np.nonzero(receptor_expression_data)[0] # Get the array of non-zero indices for receptor expressions

    # Before proceeding, we will construct the target gene matrices, as described in Wang et al. (2019), if they have been specified
    if upregulated_genes:

        num_upregulated_genes = len(upregulated_genes) # Number of upregulated target genes specified

        upregulated_gene_expressions = csr_matrix((num_cells, num_upregulated_genes)) # Initialise sparse matrix

        # Fill the matrices with the gene expressions
        for i in range(num_upregulated_genes):

            gene = upregulated_genes[i] # Get the target gene ID

            if gene in gene_ids: # If there is expression data for this gene

                gene_column_index = gene_ids.get_loc(gene)
                gene_expression_levels = gene_expression_data[:, gene_column_index]
                gene_expression_levels = csr_matrix(gene_expression_levels).transpose() # This makes sure we're assigning a vector of the correct size

                if (gene_expression_levels.max() != 0.0): # Only fill the matrix if there's a non-zero gene expression

                    upregulated_gene_expressions[:, i] = gene_expression_levels / gene_expression_levels.max() # Normalise the expression

        # Take the column means to get the mean target gene expressions for each cell
        mean_upregulated_gene_expressions = upregulated_gene_expressions.mean(1)

    # Construct the target gene matrix of downregulated genes, if specified.
    if downregulated_genes:

        num_downregulated_genes = len(downregulated_genes) # Number of downregulated target genes

        downregulated_gene_expressions = csr_matrix((num_cells, num_downregulated_genes))

        # Fill the matrices with the gene expressions
        for i in range(num_downregulated_genes):

            gene = downregulated_genes[i] # Get the target gene ID

            if gene in gene_ids: # If there is expression data for this gene

                gene_column_index = gene_ids.get_loc(gene)
                gene_expression_levels = gene_expression_data[:, gene_column_index]
                gene_expression_levels = csr_matrix(gene_expression_levels).transpose() # This makes sure we're assigning a vector of the correct size

                if (gene_expression_levels.max() != 0.0): # Only fill the matrix if there's a non-zero gene expression

                    downregulated_gene_expressions[:, i] = gene_expression_levels / gene_expression_levels.max() # Normalise the expression

        # Take the column means to get the average target gene expression for each cell
        mean_downregulated_gene_expressions = downregulated_gene_expressions.mean(1)

    # If both the ligand and receptor expressions have a positive maximum
    if ( (np.max(ligand_expression_data) != 0.0) & (np.max(receptor_expression_data != 0.0)) ):
        ligand_expression_data /= np.max(ligand_expression_data) # Normalise so the maximum ligand expression level is 1.0
        receptor_expression_data /= np.max(receptor_expression_data) # Normalise so the maximum receptor expression level is 1.0

        # Based on the method from Wang et al. (2019), we need only consider the indices where both the ligand and receptor expressions are non-zero.
        # Let us go through the non-zero ligand and receptor indices, until we find a better method.
        for ligand_index in nonzero_ligand_indices: # This corresponds to the row index

            nonzero_column_indices = [] # To maintain sparsity, we'll only normalise at the places of non-zero probabilities
            for receptor_index in nonzero_receptor_indices: # This corresponds to the column indices

                ligand_expression = ligand_expression_data[ligand_index] # The second index needs to be specified because the expression data's an n_cell x 1 matrix (weird, I know)
                receptor_expression = receptor_expression_data[receptor_index] # Similar deal as above, but for the receptor expressions

                if (ligand_expression != 0.0)&(receptor_expression != 0.0): # We need both to be non-zero to have a non-zero probability
                    
                    alpha = np.exp(- 1.0 / (ligand_expression * receptor_expression)) # Calculate the ligand-receptor interaction probability

                    # If any upregulated target genes have been specified, we need to adjust the probability based on their expression levels
                    if upregulated_genes: # If not empty

                        # Get the mean upregated gene expression for the receptor cell
                        mean_upregulated_gene_expression = mean_upregulated_gene_expressions[receptor_index, 0] # Index this way as the mean values is an n x 1 matrix

                        # We take the exponential of the average gene expression
                        # and calculate the weighting coefficient that adjusts the probability according to
                        # the average upregulated gene expression level

                        if (mean_upregulated_gene_expression != 0.0):

                            beta = np.exp( - 1.0 / mean_upregulated_gene_expression) # Quantifies the average expression of upregulated downstream genes
                            Kappa = alpha / (alpha + beta) # Adjusts the original ligand-receptor probability by the average expression

                        else:

                            beta = 0.0
                            Kappa = 0.0

                    else: # If we don't have any information on the upregulated genes, we don't adjust the communication probability
                        
                        beta = 1.0 
                        Kappa = 1.0

                    Kappa = alpha / (alpha + beta) # Calculate the weighting coefficient that adjusts the ligand-receptor interaction by the upregulated gene expressions

                    # If any downregulated target genes have been specified, we need to adjust the probability based on their expression levels
                    if downregulated_genes: #If not empty

                        mean_downregulated_gene_expression = mean_downregulated_gene_expressions[receptor_index, 0] # Get the mean gene expression for the receptor cell

                        # We take the exponential of the average gene expression
                        # and calculate the weighting coefficient that adjusts the probability according to
                        # the average downregulated gene expression level
                        gamma = np.exp( - mean_downregulated_gene_expression) # Quantifies the average expression of downregulated downstream genes
                        Lambda = alpha / (alpha + gamma) # Adjusts the original ligand-receptor probability by the average expression

                    else: # If we don't have any information on the downregulated genes, we don't adjust the communication probability
                        
                        gamma = 1.0 
                        Lambda = 1.0

                    # Calculate the cell-cell probability (un-normalised), as per Wang et al. (2019)
                    calculated_probability = alpha * beta * gamma * Kappa * Lambda 

                    # If the probability is non-zero, we add this to the list of column indices we'll need to normalise over
                    if (calculated_probability != 0.0):

                        nonzero_column_indices.append(receptor_index)
                        communication_probabilities[ligand_index, receptor_index] = calculated_probability # Update the probability matrix

            if (communication_probabilities[ligand_index, :].sum() != 0.0): # Normalise the row if we have non-zero probabilities

                for receptor_index in nonzero_receptor_indices: # Only normalise at the non-zero entries

                    communication_probabilities[ligand_index, receptor_index] /= communication_probabilities[ligand_index, :].sum() # Normalise the probabilities so that they sum across the row to be one.

    # Normalise the matrix to make sure the probabilities are all a max of 1.
    communication_probabilities /= communication_probabilities.max()

    # As this was done in the original code, we truncate all normalised probabilities below 1e-6, as suggested by Shuxiong to filter out low-level interactions
    nonzero_indices = communication_probabilities.nonzero()
    nonzero_row_indices = nonzero_indices[0]
    nonzero_column_indices = nonzero_indices[1]

    for i in range(len(nonzero_row_indices)):
        row_index = nonzero_row_indices[i]
        column_index = nonzero_column_indices[i]

        if (communication_probabilities[row_index, column_index] < 1e-6): # Truncate if the probability is too small. This factor is chosen to match the sparsity obtained in Wang et al. (2019)
            communication_probabilities[row_index, column_index] = 0.0

            # Should probably re-normalise the communication probabilities!
            communication_probabilities[row_index, :] /= communication_probabilities[row_index, :].sum() # Makes sure they sum to one.

    return communication_probabilities # Return the matrix

def calculate_communication_probabilities_list(ann_data, ligand_receptor_pairs, upregulated_genes = [], downregulated_genes = []):
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
    Returns
    ------
    communication_probabilities (list)
        A list of sparse n_cell x n_cell communication matrices P, where P_{ij} is the probability that cell i communicates
        with cell j for a ligand-receptor pair.
    """

    # Calculate the communication probabilities for each ligand-receptor pair
    individual_probabilities = [calculate_communication_probabilities(ann_data, pair, upregulated_genes, downregulated_genes) for pair in ligand_receptor_pairs]
        
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

    for row in nonzero_rows: # Sweep by rows
        normalised_aggregated_probability_matrix[row, :] /= normalised_aggregated_probability_matrix[row, :].sum() # Normalise and we done here

    return normalised_aggregated_probability_matrix