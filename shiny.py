#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# a chopped up and sewn back together, frankensteinian version of jacobs heatmap code

"""

THINGS TO DO FOR SHINY 

- comma seperted list of subnetworks to include
- option for significance floor AND ceiling 
- S/A summary vs all terms
- T/F for clustered or not 
- figure out gene hover over (that might just be in R tho)


start with just SUMMARY VS ALL and NLOGP FLOOR AND CEILING 
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

subnetworks = ["C0_U", "C0_D", "C1_U", "C1_D", "C2_U", "C2_D", "C3_U", "C3_D", 
               "C4_U", "C4_D", "C5_U", "C5_D", 
               "P0_U", "P0_D", "P1_U", "P1_D", "P2_U", "P2_D", "P3_U", "P3_D", "P4_U", "P4_D", 
               "S0_U", "S0_D","S1_U", "S1_D", "S2_U", "S2_D", "S3_U", "S3_D", "S4_U", "S4_D"]

def clustered_data(): 

    df_summary = pd.DataFrame()
    for subnetwork in subnetworks: 
        data_path = f"/Users/laurenbell/Desktop/heatmap_shiny/data/{subnetwork}/metascape_result.xlsx"
        xlsx = pd.ExcelFile(data_path)
        
        # JUSTIFICATION FOR THE TRY/EXCEPT LOOP: try to check if theres an enrichment file, if ntot, note it for now
        # post double checking: it worked!! 
        try:
            current = pd.read_excel(xlsx, sheet_name="Enrichment")
            #print("IT WORKED: " + subnetwork)
            current["Subnetwork"] = subnetwork
            
            df_summary = pd.concat([df_summary, current[current["GroupID"].str.contains("Summary")]])
        except: 
            #print("No enriched terms: " + subnetwork + " :(")
            continue
        
    #print("you should have clustered data now!")
    return df_summary

def unclustered_data():
    df_summary = pd.DataFrame()
    for subnetwork in subnetworks: 
        data_path = f"/Users/laurenbell/Desktop/heatmaps/shiny/data/{subnetwork}/metascape_result.xlsx"
        xlsx = pd.ExcelFile(data_path)
        
        try:
            current = pd.read_excel(xlsx, sheet_name="Enrichment")
            #print("IT WORKED: " + subnetwork)
            current["Subnetwork"] = subnetwork
            
            df_summary = pd.concat([df_summary, current[current["GroupID"].str.contains("Member")]])
        except: 
            #print("No enriched terms: " + subnetwork + " :(")
            continue
        
    #print("you should have unclustered data now!")
    return df_summary


# =============================================================================
# CLUSTERING FUNCTIONS WITH DENDROGRAMS
# =============================================================================
def create_clustered_bubble_plot(data, save_path=None, figsize=(16, 8), clustered=True,
                                font_size=8, title="GO Term Enrichment with Dendrograms"):
    
    #makes a heatmap w/ bubbles and dendrograms clustering both axes

    # pivot table for both axes (from long to wide format)
    pivot_nlogp = data.pivot_table(index='Subnetwork', columns='Description', 
                                   values='nLogP', fill_value=0)
    #pivot_nlog10.to_csv('/Users/laurenbell/Desktop/' + title + '_nlog10.csv')
    pivot_genes = data.pivot_table(index='Subnetwork', columns='Description', 
                                  values='NumGenesInGOList', fill_value=0)
    #pivot_genes.to_csv('/Users/laurenbell/Desktop/' + title + '_genes.csv')

    # hierarchical clustering for the rows (gene lists) and columns (GO terms)
    # pdist = distance calculation matrix 
    # linkage = linkage matrix using euclidean distance (regular distance)
    # leaves_list = returns list of nodes after clustering for labels
    if len(pivot_nlogp.index) > 1:
        row_linkage = linkage(pdist(pivot_nlogp.values, metric='euclidean'), method ="ward") 
        row_order = pivot_nlogp.index[leaves_list(row_linkage)]
        
        rlink = pd.DataFrame(row_linkage, columns=['row1', 'row2', 'distance', 'num'])
        rlink.to_csv("row_linkage.csv")

        #row_linkage = rlink
    else:
        row_linkage = None
        row_order = pivot_nlogp.index
    
    if len(pivot_nlogp.columns) > 1:
        col_linkage = linkage(pdist(pivot_nlogp.values.T, metric='euclidean'),method="ward")
        col_order = pivot_nlogp.columns[leaves_list(col_linkage)]
        
        #clink = pd.DataFrame(col_linkage, columns=['col1', 'col2', 'distance', 'num'])
        #clink.to_csv('/Users/laurenbell/Desktop/col_linkage.csv')
    else:
        col_linkage = None
        col_order = pivot_nlogp.columns
    
    # sort both to make sure they show up cleanly
    pivot_nlogp_ordered = pivot_nlogp.loc[row_order, col_order]
    #pivot_nlog10_ordered.to_csv(title + "_nlog10q.csv")
    pivot_genes_ordered = pivot_genes.loc[row_order, col_order]
    #pivot_genes_ordered.to_csv(title + "_genes.csv")
    
    # make figure with main heatmap, two dendrograms, and the legends
    fig = plt.figure(figsize=figsize)
    
    # have to add a manual check since different # of terms need different sizes and ratios
    if clustered==True:
        gs = fig.add_gridspec(2, 2, width_ratios=[0.05, 1], height_ratios=[0.25, 1], hspace=0.01, wspace=0.055)
    else:
        gs = fig.add_gridspec(2, 2, width_ratios=[0.05, 1], height_ratios=[0.25, 1], hspace=0.01, wspace=0.04)


    ax_main = fig.add_subplot(gs[1, 1])
    
    # lists of scatterplot detail for the data
    y_pos = []
    x_pos = []
    colors = []
    sizes = []
    
    for i, gene_list in enumerate(pivot_nlogp_ordered.index):
        for j, description in enumerate(pivot_nlogp_ordered.columns):
            if pivot_nlogp_ordered.iloc[i, j] > 0:
                y_pos.append(i)
                x_pos.append(j)
                colors.append(pivot_nlogp_ordered.iloc[i, j])
                size = 20 + (pivot_genes_ordered.iloc[i, j] / pivot_genes_ordered.values.max()) * 180
                sizes.append(size)
    
    #plotting scatter data w/ (GO term, subnetwork) as plotting coords
    scatter = ax_main.scatter(x_pos, y_pos, c=colors, s=sizes, 
                             cmap="crest", alpha=0.8, edgecolors='white', linewidth=0.5)
    
    #gridlines and labels
    ax_main.set_xticks(range(len(pivot_nlogp_ordered.columns)))
    ax_main.set_yticks(range(len(pivot_nlogp_ordered.index)))
    ax_main.set_xticklabels(pivot_nlogp_ordered.columns, rotation=45, ha='right', fontsize=font_size)
    ax_main.set_yticklabels(pivot_nlogp_ordered.index, fontsize=font_size+6)
    ax_main.set_xlim(-0.5, len(pivot_nlogp_ordered.columns) - 0.5)
    ax_main.set_ylim(-0.5, len(pivot_nlogp_ordered.index) - 0.5)
    ax_main.grid(True, alpha=0.3)
    ax_main.set_axisbelow(True)
    
    # adding dendrograms on both axes 
    if row_linkage is not None:
        ax_row = fig.add_subplot(gs[1, 0])
        from scipy.cluster.hierarchy import dendrogram
        #flatr = fcluster(row_linkage, t=0.5, criterion='distance')
        dendrogram(row_linkage, distance_sort = False, orientation='left', ax=ax_row, 
                  color_threshold=0, above_threshold_color='black')
        ax_row.set_xticks([])
        ax_row.set_yticks([])
        for spine in ax_row.spines.values():
            spine.set_visible(False)
    
    if col_linkage is not None:
        ax_col = fig.add_subplot(gs[0, 1])
        #flatc = fcluster(col_linkage, t=0.5, criterion='distance')
        dendrogram(col_linkage, distance_sort = False, orientation='top', ax=ax_col,
                  color_threshold=0, above_threshold_color='black')
        ax_col.set_xticks([])
        ax_col.set_yticks([])
        for spine in ax_col.spines.values():
            spine.set_visible(False)
    
    nlog10_values = np.array(colors)
    #color_bins = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]
    #color_labels = ["0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2", "2.25", "2.5"]

    color_bins = [1.3, 1.5, 2, 2.5, 3, 3.5, 4]
    color_labels = ["1.3", "1.5", "2", "2.5", "3", "3.5", "4"]
    
    #coloring/normalizing legend
    cmap = sns.color_palette("crest", as_cmap=True)
    norm = plt.Normalize(vmin=nlog10_values.min(), vmax=nlog10_values.max())
    
    color_handles = []
    for val, label in zip(color_bins, color_labels):
        color = cmap(norm(val))
        handle = plt.scatter([], [], s=100, c=[color], alpha=0.8, 
                           edgecolors='white', linewidth=0.5)
        color_handles.append(handle)
    
    
    max_genes = pivot_genes_ordered.values.max()
    size_values = [2, 6, 10, 12, 14, 18] 
    size_labels = ["2", '6','10', '14', "18"]
    size_handles = []
    
    for size_val, label in zip(size_values, size_labels):
        scaled_size = 20 + (size_val / max_genes) * 180
        handle = plt.scatter([], [], s=scaled_size, c='gray', alpha=0.6, 
                           edgecolors='white', linewidth=0.5)
        size_handles.append(handle)
    
    # legends - size and color
    if clustered==True:
        c_location=(1.01, .85)
    else:
        c_location=(1.01, .85)

    color_legend = ax_main.legend(color_handles, color_labels, bbox_to_anchor=c_location,
                                title='-LogP Value', loc='upper left', frameon=True,
                                fontsize=font_size+4, title_fontsize=font_size+3)
    color_legend.get_frame().set_facecolor('white')
    color_legend.get_frame() #.set_alpha(0.8)

    if clustered==True:
        s_location=(1.01, .35)
    else:
        s_location=(1.01, .45)

    size_legend = ax_main.legend(size_handles, size_labels, bbox_to_anchor=s_location,
                                title='Number of Genes', loc='upper left', frameon=True,
                                fontsize=font_size+4, title_fontsize=font_size+3)
    size_legend.get_frame().set_facecolor('white')
    size_legend.get_frame() #.set_alpha(0.8)

    ax_main.add_artist(color_legend)
   
    ax_main.text(1.01, .9, 'Subnetwork Labels:\nU = Upregulated\nD = Downregulated', 
                 transform=ax_main.transAxes,
                 fontsize=font_size+2, alpha=0.8, weight='normal')

    fig.suptitle(title, fontsize=font_size+15, y=0.95)
    
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
    
    return fig

# =============================================================================
# ONLY SHOWS THE TOP n TERMS IN EACH LIST
# =============================================================================
def filter_top_GO_terms_by_GeneList(data, top_n=4, nlogp_c=1.0, nlogp_f=5.0):
    """
    Filters the data for the top N GO terms from each GeneList based on nLog10_q values.
    
    Parameters:
    - data: DataFrame containing the data.
    - top_n: Number of top GO terms to keep for each GeneList, based on nLog10_q values.
    
    Returns:
    - DataFrame with filtered data.
    """
    
    data["nLogP"] = data["LogP"] * -1
    data = data[(data["nLogP"] >= nlogp_f) & (data["nLogP"] <= nlogp_c)]

    data["NumGenesInGOList"] = data["Symbols"].str.split(',').str.len()
    #print(data["NumGenesInGOList"])

    data["Subnetwork"] = data["Subnetwork"].str.replace('_', ' - ', regex=False)

    # Sort the data by GeneList and nlog10_q in descending order to prioritize higher nlog10_q values
    sorted_data = data.sort_values(by=["Subnetwork", "nLogP"], ascending=[True, False])
    
    # Group by GeneList and keep only the top_n rows per group
    top_GO_terms = sorted_data.groupby("Subnetwork").head(top_n)
    
    return top_GO_terms

def startup(clustered, nlpc, nlpf):
  if clustered == 1: 
        data = clustered_data()
        #data.to_csv("V5/clustered_data.csv", index=False)
        
        filtered_data = filter_top_GO_terms_by_GeneList(data, top_n=15, nlogp_c=nlpc, nlogp_f=nlpf)
        #filtered_data.to_csv("filtered_clustered_data.csv")

        create_clustered_bubble_plot(
            filtered_data, 
            save_path=f'/Users/laurenbell/Desktop/heatmap_shiny/top_15_clustered.png',
            figsize=(35,8),
            clustered=True,
            font_size=9,
            title=f"Top 15 Summary GO Terms per Subnetwork and Regulation"
        )

  if clustered == 2: 
      data = unclustered_data()
      #data.to_csv("V5/unclustered_data.csv", index=False)

      filtered_data = filter_top_GO_terms_by_GeneList(data, top_n=15, nlogp_c=nlpc, nlogp_f=nlpf)
      #filtered_data.to_csv("filtered_unclustered_data.csv")

      create_clustered_bubble_plot(
            filtered_data, 
            save_path=f'/Users/laurenbell/Desktop/heatmaps/shiny/top_15_unclustered.png',
            figsize=(55,10),
            clustered=False,
            font_size=9,
            title=f"Top 15 GO Terms per Subnetwork and Regulation"
      )

def main():
    clustered = 1
    nlpc=5
    nlpf=1
    startup(clustered, nlpc, nlpf)

if __name__ == "__main__":
    main()
