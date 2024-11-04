# Functions
import re
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# general and small
def space():
    print("\n")
    print("----------------")
    print("\n")


# 2_h5ad_input.py

def filter_adata_by_criteria(adata, criteria):
    """
    Filters the AnnData object based on the criteria specified in the dictionary.
    
    Parameters:
    adata (AnnData): The input AnnData object to be filtered.
    criteria (dict): A dictionary with keys 'donor', 'technique', and 'anatomy', each containing a list of strings to filter on.
    
    Returns:
    np.ndarray: An array of unique values in the 'manip' column after applying the filters.
    """
    donor_criteria = criteria['donor']
    technique_criteria = criteria['technique']
    anatomy_criteria = criteria['anatomy']
    
    # Construct the filter mask
    donor_mask = adata.obs['manip'].str.contains('|'.join(donor_criteria))
    technique_mask = adata.obs['manip'].str.contains('|'.join(technique_criteria))
    anatomy_mask = adata.obs['manip'].str.contains('|'.join(anatomy_criteria))
    
    combined_mask = donor_mask & technique_mask & anatomy_mask
    
    # Filter the adata object and return unique values from the 'manip' column
    filtered_manip = np.unique(adata.obs[combined_mask]['manip'])
    
    # return list
    return filtered_manip

## To get proper barcodes
def split_cell_id(s):
    pattern = r'([A-Z0-9]+_[A-Z0-9]+_(NAS|PRO|INT|NAS1|PRO1|INT1|TRAC|TRAC1))([ATGC]+)([0-9-]+)'
    match = re.match(pattern, s)
    if match:
        return [match.group(1), match.group(2), match.group(3)]
    else:
        return None

# selecting only top cells
def subset_by_decile(celltype_adata, celltype, manip, plots_output_dir):

    sorted_indices = celltype_adata.obs['total_counts'].sort_values(ascending=False).index

    # Calculate decile positions
    num_cells = len(celltype_adata.obs)

    check = 0
    
    if(num_cells >= 100):
        # calculate the dimension of a decile
        # and by num_cells - dim_decile we know 
        # if after selecting from 2nd decile we have
        # still at least 100 cells
        dim_decile = int(num_cells/10)

        # we need to check if the decile is > 100,
        # otherwise we will have < 100 for no reason.
        # for example:
        # if num_cells = 294, decile_1 = 29, decile_2 = 58
        # which means that we will have only 29 cells
        if dim_decile >= 100 :
            # we can normally calculate the decile
            decile_1 = int(np.floor(0.1 * num_cells))
            decile_2 = int(np.floor(0.2 * num_cells))
            # Ensure to keep only the top 100 cells from the second decile
            second_decile_indices = sorted_indices[decile_1:decile_2]
            top_100_second_decile_indices = second_decile_indices[:100]
            
        else :
            # If we do not have at least 100 cells in a decile
            # for example if we had 200 cells, 
            # each decile would contain 20 cells
            # 200 - 20 is 180. 
            # We can start taking cells from the 2nd decile
            # from index 19 to index 119 (excluded)
            # and still have 100 cells in our results.
            # We calculate check to know if we have at least 100 cells
            # if we start taking from beginning of 2nd decile
            check = num_cells - dim_decile 
            if check >= 100:
                # Extract the cells in the second decile
                decile_1 = int(np.floor(0.1 * num_cells))
                decile_2 = decile_1 + 100
                top_100_second_decile_indices = sorted_indices[decile_1:decile_2]
                
                
            else:
                # If check < 100, it means that starting from 
                # the beginning of the 2nd decile will not assure us to
                # to have at least 100 cells. 

                # we have for example 110 cells
                # which is more than 100
                # but if 110/10 = 11
                # 110 - 11 = 99
                # meaning we will not have enough cells
                # if we use the the deciles reasoning.
                # for this reason we will simply keep cells in the exact middle.
                # In this example from 4 to 94 
                bornes = int(num_cells - 100)

                # until 102 num_cells
                # with 101 or 100 it would not work
                if bornes >= 2:
                    decile_1 = int(bornes/2 - 1)
                    decile_2 = int(decile_1 + 100)

                # if I have 100 or 101 num_cells
                else:
                    decile_1 = 0
                    decile_2 = 100
                
                print("decile_1 ", decile_1)
                print("\n")
                print("decile_2 ", decile_2)

                top_100_second_decile_indices = sorted_indices[decile_1:decile_2]

    elif num_cells >= 1:
        decile_1 = 0
        decile_2 = num_cells

        top_100_second_decile_indices = sorted_indices[decile_1:decile_2]

    
    print(check)
    print("\n")
    print("\n")
    print("\n")
    print("decile_1 ", decile_1)
    print("\n")
    print("decile_2 ", decile_2)
    print("\n")
    print("\n")
    print("\n")

    # will filter after plotting
    ##################

    # Plot UMI of these best cells compared to the others
    # Create a bar plot of UMI counts per cell
    plt.figure(figsize=(10, 5))
    plt.bar(range(num_cells), celltype_adata.obs['total_counts'].sort_values(ascending=False))
    plt.xlabel('Cells')
    plt.ylabel('UMI Counts')

    # If we took real decile
    if check >= 100:
        # Add two red lines to indicate the positions of the selected 100 cells
        plt.axvline(x=decile_1, color='r', linestyle='--', label='Start of 2nd Decile')
        plt.axvline(x=decile_1 + min(100, decile_2 - decile_1), 
                    color='r', linestyle='--', label='End of Top 100 in 2nd Decile')
    else:
        # Add two red lines to indicate the positions of the selected 100 cells
        plt.axvline(x=decile_1, color='r', linestyle='--', label='Start of selected cells')
        plt.axvline(x=decile_1 + decile_2, 
                    color='r', linestyle='--', label='End of selected cells')

    # Add title
    plt.title('adata UMI per cell before filtering for celltype ' + celltype + ' in manip ' + manip, fontsize=15)
    # Add legend
    plt.legend()
    # Save the plot to a PDF file
    plt.savefig(plots_output_dir + '/UMI_counts_per_cell_' + celltype + '_' + manip + '.pdf')
    # Show the plot
    plt.show()

    ##################
    # Filtering
    ###################
    # Filter the AnnData object to keep only the selected cells
    return top_100_second_decile_indices




## save cell IDs for each Epithelial cell type
def save_cellID(celltype, celltype_adata, csv_dir, manip, plots_output_dir):

    # we filter to keep only the 100 top cells for each celltype
    # but we want also to keep trace of that because at the 
    # end we will have a global adata with all the 100 top 
    # for all cellytpes for all manip
    top_100_second_decile_indices = subset_by_decile(celltype_adata, celltype, manip, plots_output_dir)
    print("How many top_100_second_decile_indices? ", len(top_100_second_decile_indices))
    celltype_adata = celltype_adata[top_100_second_decile_indices]
    print(celltype_adata)

    # manip es: "D495_BIOP_INT"
    # dictionary of lists
    dict_celltype = {'obs_names': celltype_adata.obs_names.to_list()}
    df_celltype = pd.DataFrame(dict_celltype)

    # RENAME TO SAVE TO CSV (NOT FOR H5AD)
    print("BEFORE FILTERING: ")
    print(celltype_adata.obs_names[1:10])
    celltype_adata.obs_names = [el[2] for el in celltype_adata.obs_names.map(split_cell_id)]
    df_celltype["obs_names"] = celltype_adata.obs_names
    df_celltype["obs_names"] = df_celltype["obs_names"] + "-1"
    df_celltype["celltype"] = celltype + "_" + manip
    print("AFTER FILTERING: ")
    print(df_celltype["obs_names"][1:10])

    # saving the dataframe

    df_celltype.to_csv(csv_dir + "/epithelial_ids/" + celltype + "_cellIDs_discovair_v11_" + manip + "_epithelial.csv", 
                       index=False, header=False)
    
    return top_100_second_decile_indices


## subset the epithelial dataset for each manip, for each celltype
def subset_epithelial_celltypes(manip, adata, lst_celltypes_epithelial, csv_dir, plots_output_dir):
    print("Starting analysis for: ", manip)
    print("\n")

    lst_manip = []

    # Analysis
    # es manip : 'D495_BIOP_INT'
    adata = adata[adata.obs['manip'] == manip]

    ## subset data
    for celltype in lst_celltypes_epithelial:
        celltype_adata = adata[adata.obs.celltype_lv1_V6 == celltype]

        # if there is no data for the celltype
        if celltype_adata.n_obs == 0:
            print("No data for celltype: ", celltype)
            print("\n")
            print("----------------")
            # continue
        
        else :
            print("Extracting barcodes for: ", celltype)
            print("\n")
            print("----------------")
            print("\n")
            print("We have celltype_adata.n_obs = ", celltype_adata.n_obs, "\n")
            top_100_second_decile_indices = save_cellID(celltype, celltype_adata, csv_dir, manip, plots_output_dir)
            print("Done with celltype: ", celltype)
            print("\n")
            print("----------------")
            print("\n")

            lst_manip.extend(top_100_second_decile_indices)
        
    print("Done with all celltypes for manip: ", manip)
    print("\n")
    print("----------------")
    print("\n")
    print("----------------")
    print("\n")
    print("\n")

    return lst_manip
