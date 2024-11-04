# 4. Coverage and allelic frequencies:
#    BAMs are still in the form of CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
#    a. Calculating coverage for each position (celltype by celltype).
#    b. Calculating allelic frequencies
#    c. Calculate score_heteroplasmy from allelic frequencies
#    d. Plotting nb_reads_per_position vs score_heteroplasmy to choose cutoff of robust positions
#       We will call this value x_coverage

# To install MitoTrace: I created a conda env r_heteroplasmy_pipeline with version R version 4.0.2 
# Then I opened R session and : 
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-58.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# install.packages("seqinr")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Rsamtools")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.4.3.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# library("devtools")
# install_github("lkmklsmn/MitoTrace")

# install.packages("ggplot2")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/permute/permute_0.9-5.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/vegan/vegan_2.6-4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# install.packages("heatmaply")




library(MitoTrace)
library(ggplot2)
library(heatmaply)
# library(FactoMineR)
# library(factoextra)


# Get the current working directory
dir <- getwd()
# Get the parent directory
parentdir <- dirname(dir)
outputdir <- paste0(parentdir, "/output/4_matrices")
dataPath=paste0(parentdir, "/data/bam/2_C_celltype")
bams <- list.files(dataPath, full.names = T, pattern = ".bam$")
fasta_loc <- paste0(parentdir, "/data/fasta/GRCH38_MT.fa")

# a. Calculating coverage for each position (celltype by celltype).
mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "chrM" )
# b. Calculating allelic frequencies
af <- calc_allele_frequency(mae_res)

# To save results
# removing the extension from the name of the file
colnames(af) <- unlist(lapply(colnames(af), function(x) strsplit(x, ".", fixed=T)[[1]][1]))
af_mat <- as.matrix(af)
write.csv(af_mat, paste0(outputdir, "/af_chrM_epithelial_CELLTYPE.csv"))
print("Allelic frequencies from the coverage of celltypes, calculated.")
print("\n")

# removing the extension from the name of the file
colnames(mae_res) <- unlist(lapply(colnames(mae_res), function(x) strsplit(x, ".", fixed=T)[[1]][1]))
mae_res_mat <- as.matrix(mae_res[["read_counts"]])
write.csv(mae_res_mat, paste0(outputdir, "/mae_res_chrM_epithelial_CELLTYPE.csv"))
print("Coverage for each position, celltype by celltype, calculated.")
print("\n")

# c. Calculate score_heteroplasmy from allelic frequencies
# For each position (each 3 rows), we calculate the % of celltypes (columns), that :
# 1. first result : have a value < 20% to know the non-haplotype mutations 
# 2. second result : have a value < 0.01, to know the possible mutations of interest 

# Function to calculate the % of columns, for a row, that 
# contain a value < threshold (for example: how many celltypes 
# have a value < 0.20 of allelic frequence, for a mutation 
# at a certain position)
calculate_percentage_below_threshold <- function(row, threshold) {
  mean(row < threshold) * 100
}

# Initialize vectors for resulting % 
# These will be the vectors that will tell us 
# the % of columns for a row that respect the threshold
percent_below_20 <- numeric(nrow(af))
percent_below_1 <- numeric(nrow(af))

# For each mitochondrial genomic position (separately for each three rows)
# we want to know the % of celltype_manip (columns) where the value
# is < 20 % (not the aplotype)
# and then 1%

# Iterate through each row to calculate the %
for (i in 1:nrow(af)) {
  # extract row values
  row <- as.numeric(af[i, ])
  # calculate % of columns in row that have a value < 20% (then < 1%)
  # at the same time, we add value to the vectors for results
  percent_below_20[i] <- calculate_percentage_below_threshold(row, 0.2)
  percent_below_1[i] <- calculate_percentage_below_threshold(row, 0.01)
}

# for each position, we store the % of celltype_manip
# having a value < 20% and < 1%
# Create df with results
results_df <- data.frame(
  Position = rownames(af),
  Percent_Below_20 = percent_below_20,
  Percent_Below_1 = percent_below_1
)
print("\n")
print("For each position, we have stored the % of celltype_manip that have a value < 20% and < 1%")


# Save results
write.csv(results_df, paste0(outputdir, "/percentages_below_thresholds.csv"), row.names = TRUE)

# d. Plotting nb_reads_per_position vs score_heteroplasmy to choose cutoff of robust positions

# Assume that mae_res is the resulting matrix from MitoTrace with positions as rows and BAM files as columns
# For demonstration purposes, let's create a mock mae_res matrix
# mae_res <- matrix(runif(300, 100, 1000), nrow=100, ncol=3)  # Remove this in real usage

# Calculate the median number of reads per position (across all BAM files)
print(head(mae_res_mat))
med_reads_per_position <- apply(mae_res_mat, 1, median)

# Add the median reads per position to the data frame
results_df$Med_Reads_Per_Position <- med_reads_per_position

# Generate the plot with ggplot2
p <- ggplot(results_df, aes(x = Med_Reads_Per_Position)) +
  # Add a line for the percentage of values below 20%
  geom_line(aes(y = Percent_Below_20, color = "% celltype_manip whose af < 20%")) +
  # Add a line for the percentage of values below 0.01%
  geom_line(aes(y = Percent_Below_1, color = "% celltype_manip whose af <1%")) +
  # Add labels and title to the plot
  labs(
    title = "Median Reads per Position vs. Heteroplasmy Scores",
    x = "Median Reads per Position",
    y = "Heteroplasmy Score (%)",
    color = "Legend"
  ) +
  # Use a minimal theme for a clean look
  theme_minimal()

# Save the plot to a file
ggsave(paste0(parentdir, "/output/plots/", "plot_median_reads_vs_heteroplasmy.pdf"), plot = p)



# Iterate over each column of the matrix
for (i in seq_len(ncol(af_mat))) {
  # Select the column
  column <- af_mat[, i]
  
  # Calculate the average every 3 rows
  average <- sapply(seq(1, length(column), by = 3), function(j) mean(column[j:min(j+2, length(column))]))
  
  # Create the dataframe for ggplot
  df <- data.frame(Position = seq_along(average), Frequency = average)
  
  # Get the column name
  col_name <- colnames(af_mat)[i]
  # Remove ".bam" if it exists at the end of the column name
  col_name <- sub("\\.bam$", "", col_name)
  
  # Create the plot
  p <- ggplot(df, aes(x = Position, y = Frequency)) +
    geom_line() +
    ylim(0, max(average) * 1.1) +  # Set ylim a bit more than the maximum value of the list
    labs(title = paste("Plot for: ", col_name),
         x = "Mitochondrial genomic position",
         y = "Allelic frequency",
         caption = "For each genomic position, we calculate average of allelic frequencies (e.g. 1G>A,1G>C,1G>T).")
  
  # Save the plot
  ggsave(paste0(parentdir, "/output/plots/heteroplasmy_per_position_CELLTYPE/", "plot_heteroplasmy_", col_name,".pdf"), plot = p)
  ggsave(paste0(parentdir, "/output/plots/heteroplasmy_per_position_CELLTYPE/", "plot_heteroplasmy_", col_name,".png"), plot = p)
}


print("\n")
print("All plots done, and step 4 is done as well.")