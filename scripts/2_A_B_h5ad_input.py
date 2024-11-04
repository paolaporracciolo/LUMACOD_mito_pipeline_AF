# FIRST PART : 
# 2. Reading a h5ad as second input, with the key of the annotation of interest. 
#    (It's up to the user to filter data before using the pipeline)
#    a. Filtering the adata, for each celltype, as follows:
#       i. Sorting of cells from the one with the highest number of UMIs, to the one with the lowest one.
#       ii. Extraction of the 100 highest cells of the 2nd decile.
#           (If there are <100 cells, all cells are kept)
#       iii. Filtering the adata object keeping only the selected cells
#    b. From this filtered adata object: extraction of the barcodes in a csv file for each celltype in the form of: 
#       CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.csv

import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import re
from functions import *
import sys

###################
# INPUT VARIABLES #
h5ad_path = "/export/data2/porracciolo/mitochondria_lineage/mito_pipeline/data/h5ad/integrated_scran_V11.h5ad"
# For the moment we work with Epithelial cells, everything in the codes
# works with epithelial cells. The variable "epithelial" is not used for the moment.
# It is just a reminder I will have to do something about it.
# NB : 
# please note we filter data finding epithelial in obs.celltype_lv0_V6
# and epithelial cell types in obs.celltype_lv1_V6
epithelial = True
# manip info
criteria = {
    'donor': ["D489", "D490","D491", "D492", "D495", "D500", "D493"],
    'technique': ["BIOP"], 
    'anatomy' : ["INT", "NAS"]
}
# INPUT VARIABLES #
###################

# PATHS
dir = os.getcwd()
parentdir = os.path.dirname(dir)

# csv intermediate dir
csv_dir = os.path.join(parentdir, 'data', 'csv')
# h5ad in intermediate dir
h5ad_inter_dir = os.path.join(parentdir, 'data', 'h5ad')

# Plots output dir
plots_output_dir = os.path.join(parentdir, 'output', 'plots')

###################
# a. Filter data
# h5ad read 
adata = sc.read_h5ad(h5ad_path)
print("Correctly uploaded adata object from h5ad given path")
space()

# Preprocessing to get adata object that contains only info
# of interest
manip_lst = filter_adata_by_criteria(adata, criteria)
print("Manips: ")
print("\n")
print(manip_lst)
space()
## Selecting only data from manips of interest
adata = adata[adata.obs['manip'].isin(manip_lst)]
print("Our donors/manips adata shape: ")
print("\n")
print(adata.shape)
space()
## Selecting epithelial cells from the beginning
adata = adata[adata.obs.celltype_lv0_V6 == "Epithelial"]
print("Epithelial shape: ")
print("\n")
print(adata.shape)
space()
## Extracting list of epithelial celltypes
lst_celltypes_epithelial = adata.obs.celltype_lv1_V6.unique().to_list()
print("List of Epithelial celltypes: ")
print("\n")
print(lst_celltypes_epithelial)
space()
space()

# Filtering based on UMI
# Sort cells by the number of UMIs in descending order
# n_counts (or total_counts) is the number of UMI counts per cell
# calculate total_counts
sc.pp.calculate_qc_metrics(adata, use_raw = False)

###################
###################
# In a., we filter looking for the right manips, epithelial cells, 
# and filtered the best cells based on their number UMI

# b. From this filtered adata object: extraction of the barcodes in a csv file 
# for each celltype in the form of: 
# CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.csv

# Please note we do all of that here :
top_100_second_decile_indices_all = []

for manip in manip_lst:
    lst_manip = subset_epithelial_celltypes(manip, adata, lst_celltypes_epithelial, csv_dir, plots_output_dir)
    top_100_second_decile_indices_all.extend(lst_manip)

# adata.obs_names.intersection(top_100_second_decile_indices_all)
# and not directly top_100_second_decile_indices_all just in case
# there are indexes that do not match
print("len(top_100_second_decile_indices_all) - len(top_100_second_decile_indices_all)= ", 
      len(top_100_second_decile_indices_all) - len(top_100_second_decile_indices_all))

adata = adata[adata.obs_names.intersection(top_100_second_decile_indices_all)]
adata.write_h5ad(h5ad_inter_dir + "/adata_epithelial_topumi.h5ad")
print("barcodes are extracted")
space()