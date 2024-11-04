#!/bin/bash

# FIRST PART : 
# 3. Stats:
#    a. Counting reads
#    b. LOREM IPSUM

#    a. Counting reads
dir="$PWD"
parentdir="$(dirname "$dir")"
DATADIR=${parentdir}/data/bam/2_C_celltype
DATADIR_CELL=${parentdir}/data/bam/2_E_cell
OUTDIR=${parentdir}/"output/3_stats"

for bamfile in ${DATADIR}/*.bam; do

    echo $bamfile

    # Count the number of reads in the BAM file
    reads=$(samtools view -c $bamfile)

    echo
    echo $reads
    echo

    # Count the number of unique cells in the BAM file
    # cells=$(samtools view $bamfile | cut -f 1 | sort | uniq | wc -l)

    # Add the results to the CSV file in the output directory
    echo "$bamfile,$reads" >> $OUTDIR/nb_reads_cells_chrM_epithelial.csv
done

for bamfile in ${DATADIR_CELL}/*.bam; do

    # Count the number of reads in the BAM file
    reads=$(samtools view -c $bamfile)

    # Count the number of unique cells in the BAM file
    # cells=$(samtools view $bamfile | cut -f 1 | sort | uniq | wc -l)

    # Add the results to the CSV file in the output directory
    echo "$bamfile,$reads" >> $OUTDIR/nb_reads_cells_BARCODES_chrM_epithelial.csv
done

echo
echo "DONE"