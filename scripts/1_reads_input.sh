#!/bin/bash

# FIRST PART : 
# 1. Taking a folder containing bam/bai files as input.
#    a. Checking if all data is indexed, if not, it indexes data.
#    b. Filtering of BAMs to keep only chrM data.

###################
# INPUT VARIABLES #
DATADIR="/export/data2/porracciolo/mitochondria_lineage/mito_pipeline/data/bam/input"
# INPUT VARIABLES #
###################

dir="$PWD"
parentdir="$(dirname "$dir")"
OUTDIR=${parentdir}/data/bam/1_chrM
## Loop through each file .bam in DATADIR
for BAMFILE in ${DATADIR}/*.bam; do
    ## Extract manip name
    MANIP=$(basename ${BAMFILE} .bam)
    # a. Checking if all data is indexed, if not, it indexes data.
    ## Check if index file already exists
    if [ ! -f ${BAMFILE}.bai ]; then
        ## index data
        samtools index ${BAMFILE}
    else
        echo "Index file for ${MANIP} already exists."
    fi
    # b. Filtering of BAMs to keep only chrM data.
    ## samtools extract chrM reads
    samtools view -b ${BAMFILE} chrM > ${OUTDIR}/${MANIP}_chrM.bam
    # c. Reindex the filtered BAM file.
    samtools index ${OUTDIR}/${MANIP}_chrM.bam
done
