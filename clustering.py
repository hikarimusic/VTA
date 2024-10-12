import sys
import os
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
from sklearn.preprocessing import StandardScaler

def generate_dendrogram_heatmap(cohort_file, group_column):
    print(f"[Reading Data] ...", end='\r')
    # Construct the path to the summarize.csv file
    cohort_dir = os.path.splitext(cohort_file)[0]
    summarize_file = os.path.join(cohort_dir, 'summarize.csv')
    
    # Check if the summarize.csv file exists
    if not os.path.exists(summarize_file):
        print(f"\n[Hierarchy Cluster] Error: {summarize_file} not found")
        sys.exit(1)
    
    # Read the summarize.csv file
    df = pd.read_csv(summarize_file)
    
    # Find the index of the START_GENE column
    start_gene_index = df.columns.get_loc('START_GENE')
    
    # Extract features (gene expression data) and metadata
    gene_data = df.iloc[:, start_gene_index + 1:]  # +1 to start from the column after START_GENE
    metadata = df.iloc[:, :start_gene_index + 1]  # Include START_GENE column in metadata
    print("[Reading Data] Complete")

    # Select top n_genes most variable genes based on original data
    print(f"[Filtering Genes] ...", end='\r')
    expressed_genes = gene_data.columns[gene_data.mean() > 1]
    filtered_gene_data = gene_data[expressed_genes]
    gene_variances = filtered_gene_data.var()
    high_var_genes = gene_variances[gene_variances > 0].index
    selected_gene_data = filtered_gene_data[high_var_genes]
    print(f"[Filtered Genes] {selected_gene_data.shape[1]}")

    # Calculate correlations for genes using Pearson correlation on original values
    print(f"[Hierarchy Cluster] ...", end='\r')
    gene_dist = pdist(selected_gene_data.T, metric='correlation')
    
    # Calculate correlations for cases using Pearson correlation on original values
    case_dist = pdist(selected_gene_data, metric='correlation')
    
    # Perform hierarchical clustering on genes
    gene_linkage = hierarchy.linkage(gene_dist, method='ward')
    
    # Perform hierarchical clustering on cases
    case_linkage = hierarchy.linkage(case_dist, method='ward')
    case_order = hierarchy.leaves_list(case_linkage)
    gene_order = hierarchy.leaves_list(gene_linkage)
    
    # Reorder the data according to the clustering
    data_ordered = selected_gene_data.iloc[case_order, gene_order]
    print("[Hierarchy Cluster] Complete")
    
    # Calculate z-scores for visualization
    print("[Create Heatmap] ...", end='\r')
    scaler = StandardScaler()
    z_scores = pd.DataFrame(scaler.fit_transform(data_ordered), 
                            columns=data_ordered.columns, 
                            index=data_ordered.index)
    
    # Set up the matplotlib figure
    fig = plt.figure(figsize=(13, 11.5))
    
    # Create a gridspec for the layout
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 10, 1.5], height_ratios=[1, 0.1, 10],
                          left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.02, hspace=0.02)
    
    # Gene dendrogram
    ax_gene_dendrogram = fig.add_subplot(gs[2, 0])
    hierarchy.dendrogram(gene_linkage, orientation='left', ax=ax_gene_dendrogram, link_color_func=lambda k: 'black')
    ax_gene_dendrogram.axis('off')
    
    # Case dendrogram
    ax_case_dendrogram = fig.add_subplot(gs[0, 1])
    hierarchy.dendrogram(case_linkage, ax=ax_case_dendrogram, link_color_func=lambda k: 'black')
    ax_case_dendrogram.axis('off')
    
    # Heatmap
    ax_heatmap = fig.add_subplot(gs[2, 1])
    vmin = np.percentile(z_scores.values, 1)
    vmax = np.percentile(z_scores.values, 99)
    sns.heatmap(z_scores.T, cmap='seismic', center=0, vmin=vmin, vmax=vmax,
                xticklabels=False, yticklabels=False, cbar=False, ax=ax_heatmap)
    
    # Group indicator
    ax_groups = fig.add_subplot(gs[1, 1], sharex=ax_heatmap)
    ax_groups.set_ylim(0, 1)
    
    # Get unique groups and assign colors
    unique_groups = metadata[group_column].unique()
    color_palette = sns.color_palette("tab10", n_colors=len(unique_groups))
    color_map = dict(zip(unique_groups, color_palette))
    
    # Plot color bars for each sample's group
    for i, sample in enumerate(data_ordered.index):
        group = metadata.loc[sample, group_column]
        ax_groups.axvspan(i, i+1, facecolor=color_map[group], alpha=1)
    
    ax_groups.set_xlim(0, len(data_ordered))
    ax_groups.axis('off')
    
    # Legend and colorbar
    ax_legend = fig.add_subplot(gs[2, 2])
   
    # Add legend for groups
    legend_elements = [Rectangle((0, 0), 0.5, 0.5, facecolor=color_map[group], label=group) for group in unique_groups]
    blank_elements = [Rectangle((0, 0), 0.5, 0.5, facecolor="white", label="") for _ in range(4)]
   
    cmap = plt.get_cmap('seismic')
    z_min = round(np.floor(vmin / 0.2))
    z_max = round(np.ceil(vmax / 0.2))
    z_values = list(range(z_max, z_min, -1))
    cbar_elements = []
    for z in z_values:
        if z % 5 == 0:
            label = f'{round(z*0.2)}'
        else:
            label = ''
        cbar_elements.append(Rectangle((0, 0), 0.5, 0.5, facecolor=cmap(z*0.1/max(abs(vmax), abs(vmin))+0.5), label=label))
    
    all_elements = legend_elements + blank_elements + cbar_elements
    legend = ax_legend.legend(handles=all_elements, loc='center', 
                              ncol=1, handlelength=1, handleheight=1, 
                              handletextpad=0.5, columnspacing=0.5, labelspacing=0.0)
    legend.get_frame().set_linewidth(0.0)
    ax_legend.axis('off')
    
    # Use a clean style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Save the plot
    output_file = os.path.join(cohort_dir, f'clustering_{group_column.replace(" ", "_")}.pdf')
    plt.savefig(output_file, format='pdf', dpi=600, bbox_inches='tight')
    output_file = os.path.join(cohort_dir, f'clustering_{group_column.replace(" ", "_")}.png')
    plt.savefig(output_file, format='png', dpi=600, bbox_inches='tight')
    plt.close()
    
    print("[Create Heatmap] Complete")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 clustering.py <cohort_file> <group_column>")
        sys.exit(1)
    
    cohort_file = sys.argv[1]
    group_column = sys.argv[2]
    
    generate_dendrogram_heatmap(cohort_file, group_column)