""" 
Functions that we use to calculate the cell-cell communication probabilities. Based off the methods
from Wang et al. (2020).
"""

import numpy as np
from scipy.sparse import csr_matrix

def calculate_communication_probabilities(data, ligand_receptor_pair, upregulated_genes = [], downregulated_genes = []):
    """ 
    Method to calculate the cell-cell communication probabilities for a given ligand-receptor pair.

    Parameters
    ---------
    data (sparse csr_matrix)
        The gene expression data, with cells as the rows and genes as the columns.
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
    A sparse n_cell x n_cell communication matrix.
    """