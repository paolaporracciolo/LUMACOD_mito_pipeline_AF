# LUMACOD_mito_pipeline_AF (1st part MT_AF in Joyo)

This is the file to launch when launching the pipeline.
The pipeline allows to calculate allelic frequencies of the mitochondrial genome
across multiple bam files, for cells of a celltype annotation.
Input: 
  - h5ad file to extract barcodes of the cells
  - key of the annotation to find in adata.obs
  - folder containing bam files

The pipeline is composed of the following parts:

conda activate heteroplasmy_pipeline
## FIRST PART : 

1. Taking a folder containing bam/bai files as input.
  - a. Checking if all data is indexed, if not, it indexes data.
  - b. Filtering of BAMs to keep only chrM data.

3. Reading a h5ad as second input, with the key of the annotation of interest. 
  (It's up to the user to filter data before using the pipeline)
  - a. Filtering the adata, for each celltype, as follows:
    - i. Sorting of cells from the one with the highest number of UMIs, to the one with the lowest one.
    - ii. Extraction of the 100 highest cells of the 2nd decile.
          (If there are <100 cells, all cells are kept)
    - iii. Filtering the adata object keeping only the selected cells
  - b. From this filtered adata object: extraction of the barcodes in a csv file for each celltype in the form of:
    - CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.csv
  - c. Extract a BAM file for: 
     each celltype, of each donor, of each manip, of each compartment, in the form of:
    - CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
    - CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bai
  - d. Extract a BAM file for:
     each cell, of each celltype, of each donor, of each manip, of each compartment, in the form of:
    - CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
    - CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bai
   
5. Stats:
  - a. Counting reads
  - b. LOREM IPSUM

conda activate r_heteroplasmy_pipeline

4. Coverage and allelic frequencies:
   BAMs are still in the form of CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
  - a. Calculating coverage for each position (celltype by celltype).
  - b. Calculating allelic frequencies
  - c. Calculate score_heteroplasmy from allelic frequencies
  - d. Plotting nb_reads_per_position vs score_heteroplasmy to choose cutoff of robust positions
      We will call this value x_coverage

## SECOND PART:
Now we have the cutoff x_coverage to filter positions.

5. Filtering the matrix of allelic frequencies in 4.b., keeping positions where all celltypes
   have a coverage higher than the x_coverage cutoff. Saving these positions in robust_pos
   
6. Recalculating coverage and allelic frequencies:
   BAMs are now in the form of CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
  - a. Calculating coverage for each position (cell by cell).
  - b. Calculating allelic frequencies
  - c. Calculate score_heteroplasmy from allelic sequencies
   
7. Keeping of this second allelic frequencies matrix only robust_pos.


