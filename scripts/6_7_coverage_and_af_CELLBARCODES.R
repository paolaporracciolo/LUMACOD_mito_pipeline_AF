# 6. Recalculating coverage and allelic frequencies:
#    BAMs are now in the form of CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
#    a. Calculating coverage for each position (cell by cell).
#    b. Calculating allelic frequencies
# 7. Keeping of this second allelic frequencies matrix only robust_pos.

library(MitoTrace)
library(ggplot2)
library(heatmaply)
#library(FactoMineR)
#library(factoextra)


# Get the current working directory
dir <- getwd()
# Get the parent directory
parentdir <- dirname(dir)
outputdir <- paste0(parentdir, "/output/4_matrices")
# this is the path that changes from script 4
dataPath=paste0(parentdir, "/data/bam/2_E_cell")
bams <- list.files(dataPath, full.names = T, pattern = ".bam$")
fasta_loc <- paste0(parentdir, "/data/fasta/GRCH38_MT.fa")

# 6. Recalculating coverage and allelic frequencies:
# a. Calculating coverage for each position (cell by cell).
mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "chrM" )
# removing the extension from the name of the file
# colnames(mae_res) <- unlist(lapply(colnames(mae_res), function(x) strsplit(x, ".", fixed=T)[[1]][1]))
mae_res <- as.matrix(mae_res)
write.csv(mae_res, paste0(outputdir, "/mae_res_chrM_epithelial_BARCODES.csv"))

# b. Calculating allelic frequencies
af <- calc_allele_frequency(mae_res)
# removing the extension from the name of the file
colnames(af) <- unlist(lapply(colnames(af), function(x) strsplit(x, ".", fixed=T)[[1]][1]))
af_mat <- as.matrix(af)
write.csv(af_mat, paste0(outputdir, "/af_chrM_epithelial_BARCODES.csv"))
print("AF MAT")
print(head(af_mat))
print("\n")
print("\n")
print("\n")

# 7. Keeping of this second allelic frequencies matrix only robust_pos.
af_chrM_epithelial_CELLTYPE_above_coverage <- read.csv(paste0(outputdir, 
                                                       "/af_chrM_epithelial_CELLTYPE_above_coverage.csv"), 
                                                       header = TRUE, row.names = 1)
print("class(af_chrM_epithelial_CELLTYPE_above_coverage)")
print(class(af_chrM_epithelial_CELLTYPE_above_coverage))
print("\n")
print("class(af)")
print(class(af_chrM_epithelial_CELLTYPE_above_coverage))
print("IS THERE AN ERROR HERE?")
print("\n")
print("\n")
print("\n")

af_chrM_epithelial_BARCODES_above_coverage <- af_mat[row.names(af_mat) %in% row.names(af_chrM_epithelial_CELLTYPE_above_coverage), ]
print("OR HERE ??")
print("\n")
print("\n")
print("\n")

write.csv(af_chrM_epithelial_BARCODES_above_coverage, paste0(outputdir, "/af_chrM_epithelial_BARCODES_above_coverage.csv"), row.names = TRUE)
print("AF MA ABOVE COVERAGE")
print(head(af_chrM_epithelial_BARCODES_above_coverage))
print("\n")
print("\n")
print("\n")

print("DONE")