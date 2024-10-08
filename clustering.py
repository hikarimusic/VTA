import sys
import os
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import seaborn as sns

def generate_dendrogram_heatmap(cohort_file, group_column, start_gene, n_genes=1000):
    print(f"[Clustering] Starting analysis for {group_column}...", flush=True)
    
    # Construct the path to the summarize.csv file
    cohort_dir = os.path.splitext(cohort_file)[0]
    summarize_file = os.path.join(cohort_dir, 'summarize.csv')
    
    # Check if the summarize.csv file exists
    if not os.path.exists(summarize_file):
        print(f"[Clustering] Error: {summarize_file} not found")
        sys.exit(1)
    
    # Read the summarize.csv file
    df = pd.read_csv(summarize_file)
    
    # Find the index of the start gene
    start_gene_index = df.columns.get_loc(start_gene)
    
    # Extract features (gene expression data) and metadata
    gene_data = df.iloc[:, start_gene_index:]
    metadata = df.iloc[:, :start_gene_index]
    
    print("[Clustering] Selecting top variable genes...", flush=True)
    # Select top n_genes most variable genes based on original data
    gene_variances = gene_data.var()
    top_genes = gene_variances.nlargest(n_genes).index
    selected_gene_data = gene_data[top_genes]
    
    print("[Clustering] Calculating z-scores...", flush=True)
    # Calculate z-scores for selected genes across all samples
    z_scores = selected_gene_data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
    
    # Handle outliers by clipping z-scores
    z_scores_clipped = z_scores.clip(-4, 4)
    
    print("[Clustering] Performing hierarchical clustering on genes...", flush=True)
    # Perform hierarchical clustering on genes
    gene_linkage = hierarchy.linkage(pdist(z_scores_clipped.T), method='average')
    
    print("[Clustering] Grouping and clustering samples...", flush=True)
    # Group samples by the specified column
    grouped = metadata.groupby(group_column)
    
    # Calculate mean z-scores for each group
    group_means = z_scores_clipped.groupby(metadata[group_column]).mean()
    
    # Cluster groups based on their mean z-scores
    group_linkage = hierarchy.linkage(pdist(group_means), method='average')
    group_order = group_means.index[hierarchy.leaves_list(group_linkage)]
    
    # Cluster samples within each group
    sample_order = []
    for group in group_order:
        group_data = z_scores_clipped.loc[grouped.groups[group]]
        sample_linkage = hierarchy.linkage(pdist(group_data), method='average')
        sample_order_within_group = hierarchy.leaves_list(sample_linkage)
        sample_order.extend(group_data.index[sample_order_within_group])
    
    # Reorder the data according to the sample grouping and clustering
    z_scores_ordered = z_scores_clipped.loc[sample_order]
    
    print("[Clustering] Generating heatmap...", flush=True)
    # Set up the matplotlib figure
    fig = plt.figure(figsize=(16, 16))
    
    # Create a gridspec for the layout
    gs = fig.add_gridspec(2, 3, width_ratios=[2, 15, 1], height_ratios=[1, 15],
                          left=0.05, right=0.98, bottom=0.05, top=0.95, wspace=0.02, hspace=0.02)
    
    # Dendrogram
    ax_dendrogram = fig.add_subplot(gs[1, 0])
    
    # Compress the dendrogram by adjusting the distance threshold
    max_d = 0.7 * max(gene_linkage[:, 2])  # This sets the cut-off to 70% of the max distance
    fancy_dendrogram(
        gene_linkage,
        truncate_mode='lastp',
        p=10,  # Show only the last p merged clusters
        leaf_rotation=0.,
        leaf_font_size=12.,
        show_contracted=True,
        ax=ax_dendrogram,
        max_d=max_d,
        orientation='left'
    )
    ax_dendrogram.axis('off')
    
    # Heatmap
    ax_heatmap = fig.add_subplot(gs[1, 1])
    sns.heatmap(z_scores_ordered.T, cmap='RdBu_r', center=0, vmin=-4, vmax=4,
                xticklabels=False, yticklabels=False, cbar=False, ax=ax_heatmap)
    
    # Colorbar
    ax_colorbar = fig.add_subplot(gs[1, 2])
    norm = plt.Normalize(vmin=-4, vmax=4)
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_colorbar, orientation='vertical', label='Z-score')
    cbar.ax.set_ylabel('Z-score', fontsize=10, fontweight='bold', rotation=270, labelpad=15)
    
    # Group labels
    ax_groups = fig.add_subplot(gs[0, 1], sharex=ax_heatmap)
    ax_groups.set_ylim(0, 1)
    group_sizes = grouped.size().loc[group_order]
    start = 0
    colors = plt.cm.Set2(np.linspace(0, 1, len(group_order)))
    for i, (group, size) in enumerate(group_sizes.items()):
        end = start + size
        ax_groups.axvspan(start, end, facecolor=colors[i], alpha=0.5)
        ax_groups.text((start + end) / 2, 0.5, group, ha='center', va='center', rotation=0,
                       fontsize=10, fontweight='bold', color='black')
        start = end
    ax_groups.set_xlim(0, len(z_scores_ordered))
    ax_groups.axis('off')
    
    # Main title
    fig.suptitle(f'Gene Expression Clustering by {group_column}', fontsize=16, y=0.98, fontweight='bold')
    
    # Use a clean style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Save the plot
    output_file = os.path.join(cohort_dir, f'Dendrogram_Heatmap_{group_column.replace(" ", "_")}.pdf')
    plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[Clustering] Complete for {group_column}", flush=True)

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = hierarchy.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 clustering.py <cohort_file> <group_column> <start_gene>")
        sys.exit(1)
    
    cohort_file = sys.argv[1]
    group_column = sys.argv[2]
    start_gene = sys.argv[3]
    
    generate_dendrogram_heatmap(cohort_file, group_column, start_gene)