""" 
Functions that we use to calculate the cell-cell communication probabilities. Based off the methods
from Wang et al. (2020).
"""

import numpy as np
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
    gene_ids = gene_expression_data.var_names() # Get the list of gene IDs

    ligand = ligand_receptor_pair[0] # Get the name of the ligand
    ligand_column_index = gene_ids.get_loc(ligand) # Find where this ligand occurs along the columns of the matrix 
    ligand_expression_data = gene_expression_data[:, ligand_column_index] # Returns n_cell x 1 sparse matrix of ligand expression levels

    receptor = ligand_receptor_pair[1] # Get the name of the receptor
    receptor_column_index = gene_ids.get_loc(receptor) # Find the column index that corresponds to this receptor
    receptor_expression_data = gene_expression_data[:, receptor_column_index] # Returns n_cell x 1 sparse matrix of receptor expression levels

    # We only proceed with the probability calculations if BOTH the ligand and receptor gene expressions are non-zero for a positive number of cells.
    # This is based on the underlying assumptions of the method from Wang et al. (2019)
    nonzero_ligand_indices = np.nonzero(ligand_expression_data)[0] # Get the array of non-zero indices for ligand expressions
    nonzero_receptor_indices = np.nonzero(receptor_expression_data)[0] # Get the array of non-zero indices for receptor expressions

    # If both the ligand and receptor expressions have a positive maximum
    if ( (np.max(ligand_expression_data) > 0.0) & (np.max(receptor_expression_data > 0.0)) ):
        ligand_expression_data /= np.max(ligand_expression_data) # Normalise so the maximum ligand expression level is 1.0
        receptor_expression_data /= np.max(receptor_expression_data) # Normalise so the maximum receptor expression level is 1.0

        # Based on the method from Wang et al. (2019), we need only consider the indices where both the ligand and receptor expressions are non-zero.
        # Let us go through the non-zero ligand and receptor indices, until we find a better method.
        for ligand_index in nonzero_ligand_indices: # This corresponds to the row index
            for receptor_index in nonzero_receptor_indices: # This corresponds to the column indices

                ligand_expression = ligand_expression_data[ligand_index, 0] # The second index needs to be specified because the expression data's an n_cell x 1 matrix (weird, I know)
                receptor_expression = receptor_expression_data[receptor_index, 0] # Similar deal as above, but for the receptor expressions

                normalisation_constant = 0.0 # Initialise the normalisation constant
                if (ligand_expression > 0.0)&(receptor_expression > 0.0): # We need both to be non-zero to have a non-zero probability

                    alpha = np.exp(- 1.0 / (ligand_expression * receptor_expression)) # Calculate the ligand-receptor interaction probability

                    # If any upregulated target genes have been specified, we need to adjust the probability based on their expression levels
                    if upregulated_genes: # If not empty

                        average_upregulated_gene_expression = 0.0 # Initialise the averaged gene expression
                        num_upregulated_genes = len(upregulated_genes) # Get the number of upregulated genes specified

                        # Calculate the average gene expression level
                        for gene in upregulated_genes:
                            # Get the column index corresponding to this gene
                            gene_column_index = gene_ids.get_loc(gene)
                            gene_expression_level = gene_expression_data[receptor_index, gene_column_index]

                            average_upregulated_gene_expression += gene_expression_level

                        average_upregulated_gene_expression  /= num_upregulated_genes # Divide by the number of upstream genes to obtain the average

                            # If there's a non-zero average, we take the exponential of the average gene expression
                            # and calculate the weighting coefficient that adjusts the probability according to
                            # the average upregulated gene expression level
                        if (average_upregulated_gene_expression > 0.0): 
                            beta = np.exp( - 1.0 / average_upregulated_gene_expression) # Quantifies the average expression of upregulated downstream genes
                            Kappa = alpha / (alpha + beta) # Adjusts the original ligand-receptor probability by the average expression
                        else: # If the average gene expression level is 0.0, then beta is 0 and Kappa is 1
                            beta = 0.0
                            Kappa = 1.0

                    else: # If we don't have any information on the upregulated genes, we don't adjust the communication probability
                        
                        beta = 1.0 
                        Kappa = 1.0

                    Kappa = alpha / (alpha + beta) # Calculate the weighting coefficient that adjusts the ligand-receptor interaction by the upregulated gene expressions

                    # If any downregulated target genes have been specified, we need to adjust the probability based on their expression levels
                    if downregulated_genes: #If not empty

                        average_downregulated_gene_expression = 0.0 # Initialise the averaged gene expression
                        num_downregulated_genes = len(downregulated_genes) # Get the number of upregulated genes specified

                        # Calculate the average gene expression level
                        for gene in downregulated_genes:
                            # Get the column index corresponding to this gene
                            gene_column_index = gene_ids.get_loc(gene)
                            gene_expression_level = gene_expression_data[receptor_index, gene_column_index]

                            average_downregulated_gene_expression += gene_expression_level

                        average_downregulated_gene_expression  /= num_downregulated_genes # Divide by the number of upstream genes to obtain the average

                        # We take the exponential of the average gene expression
                        # and calculate the weighting coefficient that adjusts the probability according to
                        # the average downregulated gene expression level
                        gamma = np.exp( - average_downregulated_gene_expression) # Quantifies the average expression of downregulated downstream genes
                        Lambda = alpha / (alpha + gamma) # Adjusts the original ligand-receptor probability by the average expression

                    else: # If we don't have any information on the downregulated genes, we don't adjust the communication probability
                        
                        gamma = 1.0 
                        Lambda = 1.0

                    # Calculate the cell-cell probability (un-normalised), as per Wang et al. (2019)
                    calculated_probability = alpha * beta * gamma * Kappa * Lambda 
                    communication_probabilities[ligand_index, receptor_index] = calculated_probability

                    # Update the normalisation constant
                    normalisation_constant += calculated_probability

            communication_probabilities[ligand_index, :] /= normalisation_constant # Normalise the probabilities so that they sum across the row to be one.

    return communication_probabilities # Return the matrix