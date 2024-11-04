#!/bin/bash

dir="$PWD"
parentdir="$(dirname "$dir")"
input_folder=${parentdir}/data/bam/2_D_cell
destination_folder=${parentdir}/data/bam/2_E_cell

OUTDIR=${parentdir}/data/bam/2_C_celltype
OUTDIR_CELL=${parentdir}/data/bam/2_D_cell
## INDEX
# Loop through each file .bam in OUTDIR
for BAMFILE in ${OUTDIR}/*.bam; do
    # Extract manip name
    MANIP=$(basename ${BAMFILE} .bam)
    # index data
    samtools index ${BAMFILE}

    # Loop through each file .bam in each folder from OUTDIR
    for BAMFILE_CELL in ${OUTDIR_CELL}/${MANIP}/*.bam; do
        # Extract manip name
        MANIP_CELL=$(basename ${BAMFILE_CELL} .bam)
        # index data
        samtools index ${BAMFILE_CELL}
    done
done

## RENAME
echo $input_folder
echo

# we want to move 
for dir in $input_folder/*; do
    echo dir $dir
    echo

    if [ -d "$dir" ]; then
        # get only the folder name without the path
        dir_name=$(basename "$dir")

        echo dir_name $dir_name

        for file in "$dir"/*; do
            # get only the file name without the path
            file_name=$(basename "$file")
            # rename the file    
            mv "$file" "$dir/${dir_name}_$file_name"

        done

        # move the renamed files to the destination folder
        mv $dir/* $destination_folder

    fi
done

# remove content from 2_D_cell
# rm -r $input_folder

echo
echo "DONE"