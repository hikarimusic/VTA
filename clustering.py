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

def generate_dendrogram_heatmap(summarize_file, n_clusters, group_columns):
    print(f"[Reading Data] ...", end='\r')
    if not os.path.exists(summarize_file):
        print(f"\n[Hierarchy Cluster] Error: {summarize_file} not found")
        sys.exit(1)
    
    df = pd.read_csv(summarize_file)
    start_gene_index = df.columns.get_loc('START_GENE')
    metadata = df.iloc[:, :start_gene_index + 1]
    gene_data = df.iloc[:, start_gene_index + 1:]
    print("[Reading Data] Complete")

    # Filter and normalize genes
    print(f"[Filtering Genes] ...", end='\r')
    expressed_genes = gene_data.columns[gene_data.mean() > 0.1]
    gene_data = gene_data[expressed_genes]
    high_var_genes = gene_data.columns[gene_data.var() > 0]
    selected_gene_data = gene_data[high_var_genes]
    target_median = selected_gene_data.median(axis=1).median(axis=0)
    scale_factors = target_median / selected_gene_data.median(axis=1)
    normalized_gene_data = selected_gene_data.multiply(scale_factors, axis=0)
    print(f"[Filtered Genes] {normalized_gene_data.shape[1]}")

    # Hierarchy clustering
    print(f"[Hierarchy Cluster] ...", end='\r')
    gene_dist = pdist(normalized_gene_data.T, metric='correlation')
    case_dist = pdist(normalized_gene_data, metric='correlation')
    
    gene_linkage = hierarchy.linkage(gene_dist, method='ward')
    case_linkage = hierarchy.linkage(case_dist, method='ward')
    case_order = hierarchy.leaves_list(case_linkage)
    gene_order = hierarchy.leaves_list(gene_linkage)

    cluster_labels = hierarchy.fcluster(case_linkage, t=n_clusters, criterion='maxclust')
    
    data_ordered = normalized_gene_data.iloc[case_order, gene_order]
    print("[Hierarchy Cluster] Complete")
    
    # Create heatmap
    print("[Create Heatmap] ...", end='\r')
    scaler = StandardScaler()
    z_scores = pd.DataFrame(scaler.fit_transform(data_ordered), 
                            columns=data_ordered.columns, 
                            index=data_ordered.index)
    
    fig = plt.figure(figsize=(12.5, 11 + 0.1 * len(group_columns)))
    gs = fig.add_gridspec(2 + len(group_columns), 3, 
                          width_ratios=[1, 10, 1.5], 
                          height_ratios=[1] + [0.1] * len(group_columns) + [10],
                          left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.02, hspace=0.02)
    
    # Gene dendrogram
    ax_gene_dendrogram = fig.add_subplot(gs[1+len(group_columns), 0])
    hierarchy.dendrogram(gene_linkage, orientation='left', ax=ax_gene_dendrogram, link_color_func=lambda k: 'black')
    ax_gene_dendrogram.axis('off')
    
    # Case dendrogram
    ax_case_dendrogram = fig.add_subplot(gs[0, 1])
    hierarchy.dendrogram(case_linkage, ax=ax_case_dendrogram, link_color_func=lambda k: 'black')
    ax_case_dendrogram.axis('off')
    
    # Heatmap
    ax_heatmap = fig.add_subplot(gs[1+len(group_columns), 1])
    vmin = np.percentile(z_scores.values, 1)
    vmax = np.percentile(z_scores.values, 99)
    sns.heatmap(z_scores.T, cmap='seismic', center=0, vmin=vmin, vmax=vmax,
                xticklabels=False, yticklabels=False, cbar=False, ax=ax_heatmap)
    
    # Group indicators
    color_preset = ["Set1", "tab10", "Dark2"]
    for i, group_column in enumerate(group_columns):
        ax_groups = fig.add_subplot(gs[1+i, 1], sharex=ax_heatmap)
        ax_groups.set_ylim(0, 1)
        
        unique_groups = metadata[group_column].unique()
        color_palette = sns.color_palette(color_preset[i%3], n_colors=len(unique_groups))
        color_map = dict(zip(unique_groups, color_palette))
        
        for j, sample in enumerate(data_ordered.index):
            group = metadata.loc[sample, group_column]
            ax_groups.axvspan(j, j+1, facecolor=color_map[group], alpha=1)
        
        ax_groups.set_xlim(0, len(data_ordered))
        ax_groups.axis('off')
    
    # Legend and colorbar
    ax_legend = fig.add_subplot(gs[1+len(group_columns), 2])
   
    legend_elements = []
    for i, group_column in enumerate(group_columns):
        unique_groups = metadata[group_column].unique()
        color_palette = sns.color_palette(color_preset[i%3], n_colors=len(unique_groups))
        color_map = dict(zip(unique_groups, color_palette))
        legend_elements.extend([Rectangle((0, 0), 0.5, 0.5, facecolor=color_map[group], label=group) for group in unique_groups])
        legend_elements.extend([Rectangle((0, 0), 0.5, 0.5, facecolor="white", label="") for _ in range(4)])
   
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
    
    all_elements = legend_elements + cbar_elements
    legend = ax_legend.legend(handles=all_elements, loc='center', 
                              ncol=1, handlelength=1, handleheight=1, 
                              handletextpad=0.5, columnspacing=0.5, labelspacing=0.0)
    legend.get_frame().set_linewidth(0.0)
    ax_legend.axis('off')
    
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Save the plot
    output_dir = os.path.dirname(summarize_file)
    output_file = os.path.join(output_dir, f'clustering_{"_".join(group_columns).replace(" ", "_")}.pdf')
    plt.savefig(output_file, format='pdf', dpi=600, bbox_inches='tight')
    output_file = os.path.join(output_dir, f'clustering_{"_".join(group_columns).replace(" ", "_")}.png')
    plt.savefig(output_file, format='png', dpi=600, bbox_inches='tight')
    plt.close()
    
    # Output clustering.csv
    cohort_file = os.path.join(output_dir, 'cohort.csv')
    cohort_df = pd.read_csv(cohort_file)
    cohort_df_reordered = cohort_df.loc[data_ordered.index]
    cohort_df_reordered['cluster'] = [f"Cluster{label}" for label in cluster_labels[case_order]]

    clustering_output = os.path.join(output_dir, 'clustering.csv')
    cohort_df_reordered.to_csv(clustering_output, index=False)
    print("[Create Heatmap] Complete")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 clustering.py <cohort/summarize.csv> <group_column1> <group_column2> ... <n_clusters>")
        sys.exit(1)
    
    summarize_file = sys.argv[1]
    n_clusters = int(sys.argv[-1])
    group_columns = sys.argv[2:-1]
    
    generate_dendrogram_heatmap(summarize_file, n_clusters, group_columns)