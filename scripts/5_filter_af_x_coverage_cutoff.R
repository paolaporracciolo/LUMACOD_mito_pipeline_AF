# SECOND PART:
# Now we have the cutoff x_coverage to filter positions.
# 5. Filtering the matrix of allelic frequencies in 4.b., keeping positions where all celltypes
#    have a coverage higher than the x_coverage cutoff. Saving these positions in robust_pos

###################
# INPUT VARIABLES #
# Set the coverage cutoff
# x_coverage <- 1000  
# This gave us 0 rows (with 100% columns)

# Testing : is there at least a row where all values of coverage are higher than 0 ?
# x_coverage <- 0  
# This gave us 3 rows

x_coverage <- 1000 

# INPUT VARIABLES #
###################

# Get the current working directory
dir <- getwd()
parentdir <- dirname(dir)
outputdir <- paste0(parentdir, "/output/4_matrices")
datadir <- paste0(parentdir, "/output/4_matrices")
mae_res <- read.csv(paste0(datadir, "/mae_res_chrM_epithelial_CELLTYPE.csv"), header = TRUE, row.names = 1)
af_mat <- read.csv(paste0(datadir, "/af_chrM_epithelial_CELLTYPE.csv"), header = TRUE, row.names = 1)



# a. Filter the coverage matrix
# Identify rows where all values are greater than x_coverage
print("MAE RES")
print(head(mae_res))
print("\n")
print("\n")
print("\n")



print("MAX")
print(max(mae_res))
print("\n")
print("\n")
print("\n")



# Define the threshold percentage
threshold <- 0.6
print("threshold")
print(threshold)
print("\n")
print("\n")
print("\n")
# Modify the function to check if at least 80% of the columns meet the condition
#rows_to_keep <- apply(mae_res, 1, function(row) {
#  sum(row > x_coverage) / length(row) >= threshold
#})

rows_to_keep <- apply(mae_res, 1, function(row) {
  any(row > 0)
})
# rows_to_keep <- apply(mae_res, 1, function(row) all(row > x_coverage))
print("ROWS TO KEEP")
print(rows_to_keep)
print("\n")
print("\n")
print("\n")



print("AF MAT")
print(head(af_mat))
print("\n")
print("\n")
print("\n")



# Filter the allele frequency matrix
af_filtered <- af_mat[rows_to_keep, ]
write.csv(af_filtered, paste0(outputdir, "/af_chrM_epithelial_CELLTYPE_above_coverage.csv"), row.names = TRUE)



print("DONE")