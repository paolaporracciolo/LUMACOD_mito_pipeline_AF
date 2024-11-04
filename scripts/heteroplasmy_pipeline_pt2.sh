#!/bin/bash

# This is the second file to launch when launching the pipeline. (It's the R part)
# The pipeline allows to calculate allelic frequencies of the mitochondrial genome
# across multiple bam files, for cells of a celltype annotation.
# Input: 
#   - h5ad file to extract barcodes of the cells
#   - key of the annotation to find in adata.obs
#   - folder containing bam files

# The pipeline is composed of the following parts:

# FIRST PART : 
# 4. Coverage and allelic frequencies:
#    BAMs are still in the form of CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
#    a. Calculating coverage for each position (celltype by celltype).
#    b. Calculating allelic frequencies
#    c. Calculate score_heteroplasmy from allelic frequencies
#    d. Plotting nb_reads_per_position vs score_heteroplasmy to choose cutoff of robust positions
#       We will call this value x_coverage

# SECOND PART:
# Now we have the cutoff x_coverage to filter positions.
# 5. Filtering the matrix of allelic frequencies in 4.b., keeping positions where all celltypes
#    have a coverage higher than the x_coverage cutoff. Saving these positions in robust_pos
# 6. Recalculating coverage and allelic frequencies:
#    BAMs are now in the form of CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
#    a. Calculating coverage for each position (cell by cell).
#    b. Calculating allelic frequencies
#    c. Calculate score_heteroplasmy from allelic sequencies
# 7. Keeping of this second allelic frequencies matrix only robust_pos.

# NB: my output h5ad should have for each manip for each celltype 100 or less cells

# NB :
# For the moment you have to check in each script of this folder the section :

###################
# INPUT VARIABLES #
#
# INPUT VARIABLES #
###################

# To check the variables to change (then the user will be able to enter all of
# those variables when launching pipeline).



# conda activate r_heteroplasmy_pipeline



# FIRST TEST
# =========> command
# STATUS : did it work ?



# 4_coverage_and_af_CELLTYPE.R
# =========> nohup Rscript 4_coverage_and_af_CELLTYPE.R > 4_coverage_and_af_CELLTYPE.log &
# STATUS : WORKED with CTRL as well

# 5_filter_af_x_coverage_cutoff.R
# =========> nohup Rscript 5_filter_af_x_coverage_cutoff.R > 5_filter_af_x_coverage_cutoff.log &
# STATUS : WORKED with CTRL as well

# 6_7_coverage_and_af_CELLBARCODES.R
# =========> nohup Rscript 6_7_coverage_and_af_CELLBARCODES.R > 6_7_coverage_and_af_CELLBARCODES.log &
# STATUS : WORKED with CTRL as well