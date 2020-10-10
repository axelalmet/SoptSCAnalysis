"""
The following lines of code concern the processing of the initial single-cell RNA-seq data, whether it is 10X data, a simple txt/csv file, or even a
scanpy object. In following with the data storage conventions, the gene expression matrix itself will be stored as an AnnData, with the relevant
n_obvs cell IDs, the n_vars gene IDs.
"""

# Import the relevant methods
from _probability import calculate_communication_probabilities_list, calculate_aggregated_probability_matrix, normalise_aggregated_probabilities # Methods to calculate the communication probabilities
from _network import construct_network_from_probability_matrix, construct_network_list_from_probability_list # Methods to construct adjacency graphs

class SoptSC:
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
    num_genes (int)
        The number of genes. While we could obtain num_genes from gene_ids, this is a more convenient
        way of storing this attribute.
    cell_ids (pandas data frame)
        The cell IDs, in case one wants to annotate the resulting commmunication matrices.
    num_cells (int)
        The number of cells. We could obtain this from cell_ids, but this is just more convenient.
    signalling_pathways (dict)
        A dictionary of each signalling pathway we would like to consider, where the values are a list
        of (ligand, receptor) tuples that we associate with the signalling pathway.
    upregulated_genes (dict)
        The associated target genes that we assume/know to be up-regulated by a signalling pathway.
        Keys are the signalling pathway name, the values are a list of the target genes.
    downregulated_genes (dict)
        The associated target genes that we assume/know to be down-regulated by a signalling pathway.
        Keys are the signalling pathway name, the values are a list of the target genes.
    individual_probability_matrices (dict)
        These are the individual probability matrices calculated for each ligand-receptor pair within each
        signalling pathway. The dictionary keys are the signalling pathway names and the values are lists
        of the probability matrices. This attribute is not initialised unless specifically called.
    aggregated_probability_matrices (dict)
        These are the aggregated probability matrices, calculated as the averaged probability matrix
        over each ligand-receptor pair for a given signalling pathway.
        The dictionary keys are the signalling pathway names and the values are the aggregated matrix.
        This attribute is not initialised unless specifically called.
    individual_adjacency_matrices (dict)
        These are the adjacency matrices that correspond to the individual probabilities calculated for each
        ligand-receptor pair within each signalling pathway. We construct these using NetworkX and take care
        to maek sure these are 
    """

    def __init__(self, ann_data, signalling_pathways = {}, ligand_receptor_pairs = [], upregulated_genes = {}, downregulated_genes = {}):
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

        self.data = ann_data # Set the gene expression data in the form of an annotated data frame
        self.gene_ids = ann_data.var_names # Set the gene IDs from the annotated data matrix
        self.cell_ids = ann_data.obs # Set the cell IDs
        self.signalling_pathways = {} # Initialise the signalling pathways dictionary
        self.upregulated_genes = {} # Initialise the upregulated genes
        self.downregulated_genes = {} # Initialise the downregulated genes

        num_pathways = len(signalling_pathways)

        # Store the ligand receptor pairs with their corresponding pathways
        if ligand_receptor_pairs: # If the specified pairs are non-empty (as we can always set them later)
            
            for i in range(num_pathways):

                pathway_name = signalling_pathways[i] # Get the signalling pathway name
                pairs = ligand_receptor_pairs[i] # Get the ligand-receptor pairs corresponding to thsi pathway

                self.signalling_pathways[pathway_name] = pairs # Set the ligand-receptor pairs
        
        else: # If none have been specified, initialise a dictionary with empty values
            
            for pathway in signalling_pathways:

                self.signalling_pathways[pathway] = [] # Just initialise an empty list

        # Store the target genes if they have been specified
        if upregulated_genes:

            for i in range(num_pathways):

                pathway_name = signalling_pathways[i] # Get the signalling pathway name
                genes = upregulated_genes[i] # Get the upstream genes associated with the patywah

                self.upregulated_genes[pathway_name] = genes # Set the genes

        else: # Else we have to initialise empty lists for all of signalling pathways

            for pathway in signalling_pathways:
                self.upregulated_genes[pathway] = [] 

        if downregulated_genes:

                pathway_name = signalling_pathways[i] # Get the signalling pathway name
                genes = downregulated_genes[i] # Get the upstream genes associated with the patywah

                self.downregulated_genes[pathway_name] = genes # Set the genes

        else: # Else we initialise empty lists

            for pathway in signalling_pathways:

                self.downregulated_genes[pathway] = []


        self.individual_probabilities = {} # Initialise the individual probability matrices
        self.aggregated_probabilities = {} # Initialise the aggregated probability matrices
        self.normalised_aggregated_probabilities = {} # Initialise the normalised aggregated probability matrices
        self.individual_networks = {} # Initialise the directed networks for the individual probability matrices
        self.aggregated_networks = {} # Initialise the directed networks for the aggregated probability matrices
        self.normalised_aggregated_networks = {} # Initialise the directed networks for the normalised aggregated probability matrices

    def set_signalling_pathways(self, signalling_pathways, ligand_receptor_pairs):
        """
        Set the signalling pathways that we want to consider for the communication probabilities. 
        
        Parameters
        ----------
        signalling_pathways (list)
            The relevant signalling pathways, which form the keys for the dictionary that we 
            store everything in.
        ligand_receptor_pairs (list of lists)
            The list of ligand-receptor pairs associated with each signalling pathway.
        """

        num_pathways = len(signalling_pathways) # Get the number of signalling apthways

        # Store each pathway name and its list of ligand-reeptor pairs
        for i in range(num_pathways): 

            pathway = signalling_pathways[i] # Get the pathway name
            ligand_receptors = ligand_receptor_pairs[i] # Get the list of ligand-receptor pairs
            self.signalling_pathways[pathway] = ligand_receptors
    
    def set_upregulated_genes(self, signalling_pathways, upregulated_genes):
        """
        Set the target genes we know to be upregulated by the signalling pathways. 

        Parameters
        ----------
        signalling_pathways (list)
            The names of the signalling ways for which we want to store the upregulated gene list.
        upregulated_genes (list)
            The upregulated genes targeted by each signalling pathway.
        """

        num_pathways = len(signalling_pathways) # Get the number of signalling apthways

        # Store each pathway name and its list of upregulated genes, IF this list is non-empty
        if upregulated_genes:

            for i in range(num_pathways): 

                pathway = signalling_pathways[i] # Get the pathway name
                genes = upregulated_genes[i] # Get the list of ligand-receptor pairs
                self.upregulated_genes[pathway] = genes

        else: # This option allows us to clear the downstream genes

            for pathway in signalling_pathways:

                self.upregulated_genes[pathway] = []

    def set_downregulated_genes(self, signalling_pathways, downregulated_genes):
        """
        Set the target genes we know to be downregulated by the signalling pathways. 

        Parameters
        ----------
        signalling_pathways (list)
            The names of the signalling ways for which we want to store the upregulated gene list.
        downregulated_genes (list)
            The downregulated genes targeted by each signalling pathway.
        """

        num_pathways = len(signalling_pathways) # Get the number of signalling apthways

        # Only store the genes if there's a non-empty list of downstream genes
        if downregulated_genes:
            for i in range(num_pathways): 

                pathway = signalling_pathways[i] # Get the pathway name
                genes = downregulated_genes[i] # Get the list of ligand-receptor pairs
                self.downregulated_genes[pathway] = genes

        else: # We allow this option to clear the downstream genes
            for pathway in signalling_pathways:
                self.downregulated_genes[pathway] = []

    def calculate_individual_probabilities(self, specified_pathways):
        """
        Calculates the individual probability matrices for each ligand-receptor pair specified
        within the considered signalling pathways. These individual probabilities are then
        stored in the soptsc attribute individual_probability_matrices

        Parameters
        ----------
        specified_pathways (list)
            The specified signalling pathway we want to consider. This should correspond
            to a key in the dictionary of signalling pathways set in the SoptSC object.
        """

        ann_data = self.data # Get the annotated data frame required for probability calculations

        for pathway in specified_pathways:

            ligand_receptor_pairs = self.signalling_pathways[pathway] # Get the list of ligand-receptor pairs
            upregulated_genes = self.upregulated_genes[pathway] # Get the list of upregulated target genes
            downregulated_genes = self.downregulated_genes[pathway] # Get the list of downregulated target genes

            probability_matrices = calculate_communication_probabilities_list(ann_data, ligand_receptor_pairs, upregulated_genes, downregulated_genes)

            # Store the list
            self.individual_probabilities[pathway] = probability_matrices

    def calculate_aggregated_probabilities(self, specified_pathways):
        """
        Calculates the aggregated probability matrices for each ligand-receptor pair within a 
        provided signalling pathway.

        Parameters
        ----------
        specified_pathways (list)
            The specified signalling pathway we want to consider. This should correspond
            to a key in the dictionary of signalling pathways set in the SoptSC object.

        """

        for pathway in specified_pathways:

            # If the individual probability matrices have been calculated, we only need to set the aggregated probabilities
            if pathway in self.individual_probabilities: 

                individual_probabilities = self.individual_probabilities[pathway] # Get the individual probability matrices (one for each ligand-receptor pair)
                self.aggregated_probabilities[pathway] = calculate_aggregated_probability_matrix(individual_probabilities) # Calculate the aggregated probability from the list
    
            else: # Else we need to calculate the individual probability matrices first

                self.calculate_individual_probabilities([pathway]) # Calculate the individual probabilities

                individual_probabilities = self.individual_probabilities[pathway] # Get the just-calculated matrices
                self.aggregated_probabilities[pathway] = calculate_aggregated_probability_matrix(individual_probabilities) # Calculate the aggregated probability from the list

    def calculate_normalised_aggregated_probabilities(self, specified_pathways):
        """
        Calculates the normalised aggregated probability matrices for each ligand-receptor pair within a 
        provided signalling pathway.

        Parameters
        ----------
        specified_pathways (list)
            The specified signalling pathway we want to consider. This should correspond
            to a key in the dictionary of signalling pathways set in the SoptSC object.

        """

        for pathway in specified_pathways:

            # If the aggregated probability has already been calculated, we just normalise it.
            if pathway in self.aggregated_probabilities: 

                aggregated_probability = self.aggregated_probabilities[pathway] # Get the individual probability matrices (one for each ligand-receptor pair)
                self.normalised_aggregated_probabilities[pathway] = normalise_aggregated_probabilities(aggregated_probability) # Calculate the aggregated probability from the list
    
            else: # Else we need to calculate the aggregated probability matrices first (and maybe the individual probability matrices too!)

                self.calculate_aggregated_probabilities([pathway]) # Calculate the aggregated probabilities

                aggregated_probability = self.aggregated_probabilities[pathway] # Get the just-calculated matrices
                self.normalised_aggregated_probabilities[pathway] = normalise_aggregated_probabilities(aggregated_probability) # Calculate the aggregated probability from the list

    def construct_individual_networks(self, specified_pathways, weight_label = 'probability'):
        """
        Constructs the corresponding directed networks for individual probability matrices
        that have been calculated. These networks are stored in a similar fashion to the probabilities;
        that is, in a dictionary where the keys are signalling pathway names and the items are the lists
        of networks.

        Parameters
        ---------
        specified_pathways (list)
            The list of signalling pathways for which we want to construct the adjacency matrices. For now,
            let's assume that we've already calculated the communication probability matrices.
        weight_label (str)
            The name of the edge weights, used when constructing the graph.
        """

        for pathway in specified_pathways:

            individual_probabilities = self.individual_probabilities[pathway] # Get the list of probability matrices
            self.individual_networks[pathway] = construct_network_list_from_probability_list(individual_probabilities, weight_label)

    def construct_aggregated_networks(self, specified_pathways, weight_label = 'probability'):
        """
        Constructs the corresponding directed networks for aggregated probability matrices
        that have been calculated. These networks are stored in a similar fashion to the probabilities;
        that is, in a dictionary where the keys are signalling pathway names and the items are the lists
        of networks.

        Parameters
        ---------
        specified_pathways (list)
            The list of signalling pathways for which we want to construct the adjacency matrices. For now,
            let's assume that we've already calculated the communication probability matrices.
        weight_label (str)
            The name of the edge weights, used when constructing the graph.
        """

        for pathway in specified_pathways:

            aggregated_probability = self.aggregated_probabilities[pathway] # Get the list of probability matrices
            self.aggregated_networks[pathway] = construct_network_from_probability_matrix(aggregated_probability, weight_label)

    def construct_normalised_aggregated_networks(self, specified_pathways, weight_label = 'probability'):
        """
        Constructs the corresponding directed networks for normalised probability matrices
        that have been calculated. These networks are stored in a similar fashion to the probabilities;
        that is, in a dictionary where the keys are signalling pathway names and the items are the lists
        of networks.

        Parameters
        ---------
        specified_pathways (list)
            The list of signalling pathways for which we want to construct the adjacency matrices. For now,
            let's assume that we've already calculated the communication probability matrices.
        weight_label (str)
            The name of the edge weights, used when constructing the graph.
        """

        for pathway in specified_pathways:

            normalised_aggregated_probability = self.normalised_aggregated_probabilities[pathway] # Get the list of probability matrices
            self.normalised_aggregated_networks[pathway] = construct_network_from_probability_matrix(normalised_aggregated_probability, weight_label)

