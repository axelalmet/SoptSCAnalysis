"""
The following lines of code concern the processing of the initial single-cell RNA-seq data, whether it is 10X data, a simple txt/csv file, or even a
scanpy object. In following with the data storage conventions, the gene expression matrix itself will be stored as an AnnData, with the relevant
n_obvs cell IDs, the n_vars gene IDs.
"""

Class SoptSC:
    """
    Creates a SoptSC class with attributes needed to calculate the cell-cell communication probabilities using
    the method described in Wang et al. (2019) (https://doi.org/10.1093/nar/gkz204).

    Attributes
    -------------
    data (scipy sparse matrix)
        The gene expression data matrix, which should be sparse. The rows correspond to cells, while the
        columns correspond to the gene IDs.
    gene_ids (pandas data frame)
        The gene IDs which are used to identify which parts of the gene expression data need
        to be subsetted when calculating the cell-cell communication probabilities.
    cell_ids (pandas data frame)
        The cell IDs, in case one wants to annotate the resulting commmunication matrices.
    signalling_pathways (dict)
        A dictionary of each signalling pathway we would like to consider, where the values are a list
        of (ligand, receptor) tuples that we associate with the signalling pathway.
    upregulated_genes (dict)
        The associated target genes that we assume/know to be up-regulated by a signalling pathway.
        Keys are the signalling pathway name, the values are a list of the target genes.
    downregulated_genes (dict)
        The associated target genes that we assume/know to be down-regulated by a signalling pathway.
        Keys are the signalling pathway name, the values are a list of the target genes.
    """
    def __init__(self, ann_data, signalling_pathways = {}, upregulated_genes = {}, downregulated_genes = {}):
        """
        Initialise the SoptSC object.
        
        Parameters
        ----------
        ann_data (annData)
            An annotated data matrix which also contains the names of the gene IDs and cell IDs.
        signalling_pathways (dict)
            Signalling pathways to consider when calculating the cell-cell communication probabilities.
        upregulated_genes (dict)
            The corresponding downstream target genes that are up-regulated by the signalling pathways.
            These are used to adjust the communication probabilities.
            Note that these genes do not necessarily have to be specified.
        downregulated_genes (dict)    
            The corresponding downstream target genes that are down-regulated by the signalling pathways.
            These are used to adjust the communication probabilities.
            Note that these genes do not necessarily have to be specified.
        """
        self.data = ann_data.X # Set the raw gene expression data
        self.gene_ids = annData.var['gene_ids'] # Set the gene IDs from the annotated data matrix
        self.cell_ids = ann_data.obs # Set the cell IDs
        self.signalling_pathways = signalling_pathways # Set the signalling pathways we want to consider
        self.upregulated_genes = upregulated_genes # Set the upregulated target genes
        self.downregulated_genes = downregulated_genes # Set the downregulated target genes

    def set_signalling_pathways(self, signalling_pathways):
        """
        Set the signalling pathways that we want to consider for the communication probabilities. 
        
        Parameters
        ----------
        signalling_pathways (dict)
        The relevant signalling pathways, where keys correspond to the signalling pathway names, e.g. Wnt,
        and the values are lists of (ligand, receptor) tuples.
        """
        self.signalling_pathways = signalling_pathways
    
    def set_downregulated_genes(self, upregulated_genes):
        """
        Set the target genes we know to be upregulated by the signalling pathways. 

        Parameters
        ----------
        upregulated_genes (dict)
        The upregulated genes targeted by each signalling pathway.
        """
        self.upregulated_genes = upregulated_genes

    def set_downregulated_genes(self, downregulated_genes):
        """
        Set the target genes we know to be downregulated by the signalling pathways. 

        Parameters
        ----------
        downregulated_genes (dict)
        The downregulated genes targeted by each signalling pathway.
        """
        self.downregulated_genes = downregulated_genes

    def calculate_probability_matrix_for_ligand_receptor_pair(self, ligand_receptor_pair):
        """
        Calculates the individual probability matrix for a given ligand receptor pair
        in a signalling pathway.

        Parameters
        ----------
        ligand_receptor_pair (tuple)
            The tuple of the ligand and receptor pair we are considering. 

        Returns 
        --------
        An n_cells x n_cells sparse matrix, where n_cells is the number of cells.
        """

    def calculate_individual_probabilities_for_signalling_pathways(self, specified_pathways):
        """
        Calculates the individual probability matrices for each ligand-receptor pair specified
        within the considered signalling pathways.

        Parameters
        ----------
        specified_pathways (list)
            The specified signalling pathway we want to consider. This should correspond
            to a key in the dictionary of signalling pathways set in the SoptSC object.

        Returns
        --------
        A dictionary where the keys are the signalling pathway names and
        the values are lists of n_cell x n_cell sparse matrices, that is,
        the individual probability matrices calculated for each ligand-receptor pair
        within the signalling pathway.
        """

    def calculate_aggregated_probabilities_for_signalling_pathways(self, specified_pathways):
        """
        Calculates the aggregated probability matrices for each ligand-receptor pair within a 
        provided signalling pathway.

        Parameters
        ----------
        specified_pathways (list)
            The specified signalling pathway we want to consider. This should correspond
            to a key in the dictionary of signalling pathways set in the SoptSC object.

        Returns
        --------
        A dictionary where the keys are the signalling pathway names and
        the values are the n_cell x n_cell sparse matrices, which, for now,
        are taken as the average of the individual sparse matrices.
        """
    

