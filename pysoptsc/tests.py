from matplotlib import pyplot as plt
import scanpy as sc
from scipy.sparse import csr_matrix
from soptsc import *
from _probability import *
import networkx as nx
import collections

# First initialise some settings for scanpy
sc.settings.verbosity = 3 # Possible values: (0) errors, (1) warnings, (2) info, (3) hints
sc.settings.set_figure_params(dpi = 80, facecolor='white')

# First load the data (we have to take the transpose, because we need the cells to be the rows and genes to be the columns)
# joostdata = sc.read_text('/Users/axelalmet/Documents/scRNASeqData/Joost2016/GSE67602_Joost_et_al_expression.txt').transpose() # Directory with the text file
joostdata = sc.read_text('/Users/axelalmet/Documents/MATLAB/SoptSC/Data/JoostData.txt').transpose() # Directory with the text file
joostdata.var_names_make_unique() # If var_names = 'gene_ids', when this step isn't necessary

sc.pp.log1p(joostdata, base = 10) # For some reason Shuxiong does this

### Test that the soptsc object initialises correctly
joost_soptsc = SoptSC(joostdata)

# Test that we can store variables correctly
pathway_names = ['Tgfb', 'Wnt', 'Bmp'] # Name of the signalling_pathways
ligand_receptor_pairs = [[('Tgfb1', 'Tgfbr1'), ('Tgfb1', 'Tgfbr2'), ('Tgfb2', 'Tgfbr1'), ('Tgfb2', 'Tgfbr2')], \
                        [('Wnt3', 'Fzd1'), ('Wnt4', 'Fzd1'), ('Wnt5a', 'Fzd1'), ('Wnt6', 'Fzd1'), ('Wnt10a', 'Fzd1')], \
                        [('Bmp1', 'Bmpr2'), ('Bmp2', 'Bmpr2'), ('Bmp4', 'Bmpr2'), ('Bmp7', 'Bmpr2')]] # Name of the ligand-receptor pairs
upregulated_genes = [['Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'], \
                        ['Ctnnb1','Lgr5','Runx2','Apc','Mmp7','Dkk1','Ccnd1'], \
                        ['Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1']] # The upregulated genes

joost_soptsc.set_signalling_pathways(pathway_names, ligand_receptor_pairs) # Set the signalling pathways and their ligand-receptor pairs

# Store upregulated/downregulated genes, without the genes first
joost_soptsc.set_upregulated_genes(pathway_names, upregulated_genes)
joost_soptsc.set_downregulated_genes(pathway_names, [])

# Calculate the probabilities now for a single pathway
joost_soptsc.calculate_individual_probabilities(pathway_names, verbose=True)

# # Calculate the aggregated probabilities for the other pathways time, I want to see if it also calculates the individual matrices properly
# joost_soptsc.calculate_aggregated_probabilities(['Tgfb', 'Wnt', 'Bmp'])

# # Calculate the individual network list for a signalling pathway
# joost_soptsc.construct_individual_networks(['Tgfb'])

# # Calculated the aggregated network for a signalling pathway
# joost_soptsc.construct_aggregated_networks(['Wnt'])

