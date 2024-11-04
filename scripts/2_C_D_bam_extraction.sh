#!/bin/bash

# FIRST PART : 
# 2. Reading a h5ad as second input, with the key of the annotation of interest. 
#    c. Extract a BAM file for: 
#       each celltype, of each donor, of each manip, of each compartment, in the form of:
#       CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
#       CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bai
#    d. Extract a BAM file for:
#       each cell, of each celltype, of each donor, of each manip, of each compartment, in the form of:
#       CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
#       CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bai

###################
# INPUT VARIABLES #
# check right path in joyo for sincelor
JAVA="/export/apps/sicelore-2.1/Jar/Sicelore-2.1-jar-with-dependencies.jar" 
MANIPS=("D489_BIOP_INT" "D489_BIOP_NAS" "D490_BIOP_INT" "D490_BIOP_NAS" "D491_BIOP_INT" "D491_BIOP_NAS" "D492_BIOP_INT" "D492_BIOP_NAS" "D495_BIOP_INT" "D495_BIOP_NAS" "D500_BIOP_INT" "D500_BIOP_NAS" "D493_BIOP_INT" "D493_BIOP_NAS")
CELLTYPES=("Basal" "Suprabasal" "Cycling Basal" "Multiciliated" "Goblet" "Submucosal Glands" "Secretory" "Myoepithelial" "Squamous" "Precursor" "Deuterosomal" "Ionocyte" "PNEC")
NUMBERS=("7754" "7755" "7759" "7760" "7763" "7764" "7767" "7768" "7782" "7783" "7810" "7811" "7772" "7773")
# INPUT VARIABLES #
###################
dir="$PWD"
parentdir="$(dirname "$dir")"
DATADIR=${parentdir}/data/bam/1_chrM
CSVDIR=${parentdir}/data/csv/epithelial_ids
CSVDIR_CELL=${parentdir}/data/csv/epithelial_ids_CELL
OUTDIR=${parentdir}/data/bam/2_C_celltype
OUTDIR_CELL=${parentdir}/data/bam/2_D_cell

# Modifying names where you have empty spaces with "_"
for file in "$DATADIR"/*; do
    # Ottieni il nome base del file
    base=$(basename "$file")
    # Sostituisci gli spazi nel nome del file con underscores
    newname=$(echo $base | tr ' ' '_')
    # Rinomina il file se il nuovo nome è diverso dal vecchio
    if [ "$base" != "$newname" ]; then
        mv "$DATADIR/$base" "$DATADIR/$newname"
    fi
done

# Modifying names where you have empty spaces with "_"
for file in "$CSVDIR"/*; do
    # Ottieni il nome base del file
    base=$(basename "$file")
    # Sostituisci gli spazi nel nome del file con underscores
    newname=$(echo $base | tr ' ' '_')
    # Rinomina il file se il nuovo nome è diverso dal vecchio
    if [ "$base" != "$newname" ]; then
        echo "echo $base != $newname"
        mv "$CSVDIR/$base" "$CSVDIR/$newname"
    fi
done

# Copy all CSV files from the source directory to the destination directory
cp "$CSVDIR"/*.csv "$CSVDIR_CELL"

# Remove the second column and any trailing comma from each CSV file in the destination directory
for file in "$CSVDIR_CELL"/*.csv; do
    awk -F, '{OFS=","; $2=""; sub(/,$/, ""); print}' "$file" > temp && mv temp "$file"
done



CELLTYPES=("Basal" "Suprabasal" "Cycling_Basal" "Multiciliated" "Goblet" "Submucosal_Glands" "Secretory" "Myoepithelial" "Squamous" "Precursor" "Deuterosomal" "Ionocyte" "PNEC")

# Loop over the array MANIPS
for ((i = 0; i < ${#MANIPS[@]}; i++)); do
    # Assign the i-th element of MANIPS to MANIP/f/3512455
    MANIP=${MANIPS[$i]}
    # Assign the i-th element of NUMBERS to NB
    NB=${NUMBERS[$i]}
    # Loop over the array CELLTYPES
    for SUBSET in "${CELLTYPES[@]}"; do
        # Construct the BAM file path
        BAM="${DATADIR}/${NB}_${MANIP}_possorted_genome_chrM.bam"
        # Construct the CSV file path
        CSV="${CSVDIR}/${SUBSET}_cellIDs_discovair_v11_${MANIP}_epithelial.csv"
        CSV_CELL="${CSVDIR_CELL}/${SUBSET}_cellIDs_discovair_v11_${MANIP}_epithelial.csv"
        # Check if both BAM and CSV files exist and are not empty
        if [[ -s $BAM && -s $CSV ]]; then
            
            # C.
            # CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
            # CELLTYPE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bai
            # Construct the command to be executed
            echo
            echo
            SCRIPT="java -jar ${JAVA} SplitBamPerCluster -CSV ${CSV} -I ${BAM} -O ${OUTDIR} -CELLTAG CB"
            echo $SCRIPT
            # Execute the command
            $SCRIPT
            # Echo a message indicating the completion of this subset
            echo "Done extracting celltype bam for $SUBSET"
            echo
            echo

            # D.
            # CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bam
            # CELLTYPE_BARCODE_cellIDs_ADATANAME_DONOR_MANIP_COMPARTMENT.bai
            echo "Now doing it on all of its barcodes"
            # Command to execute the bam split not only for celltype 
            # but also for cell barcodes
            echo
            echo
            OUTDIR_CELL_celltype="${OUTDIR_CELL}/${SUBSET}_${MANIP}"
            mkdir $OUTDIR_CELL_celltype
            SCRIPT2="java -jar ${JAVA} SplitBamPerCell -CSV ${CSV_CELL} -I ${BAM} -O ${OUTDIR_CELL_celltype} -CELLTAG CB"
            echo $SCRIPT2
            # Execute the command
            $SCRIPT2
            echo "Done extracting barcode of celltype bams for $SUBSET"
            echo
            echo
        else
            # If the BAM file does not exist or is empty, echo a message
            if [[ ! -s $BAM ]]; then
                echo "The BAM file $BAM is empty or does not exist"
            fi
            # If the CSV file does not exist or is empty, echo a message
            if [[ ! -s $CSV ]]; then
                echo "The CSV file $CSV is empty or does not exist"
            fi
        fi
        # Echo an empty line for readability
        echo
        echo "______________________"
        echo "______________________"
        echo
    done
done