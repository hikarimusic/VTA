import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Rectangle

def main(summarize_file, group_column, group1, group2):
    # Read data
    print(f"[Reading Data] ...", end='\r')
    if not os.path.exists(summarize_file):
        print(f"\n[Reading Data] Error: {summarize_file} not found")
        sys.exit(1)
    
    df = pd.read_csv(summarize_file)
    start_gene_index = df.columns.get_loc('START_GENE')
    metadata = df.iloc[:, :start_gene_index + 1]
    gene_data = df.iloc[:, start_gene_index + 1:]
    print("[Reading Data] Complete")

    # Filter and normalize gene
    print(f"[Filtering Genes] ...", end='\r')
    expressed_genes = gene_data.columns[gene_data.mean() > 1]
    gene_data = gene_data[expressed_genes]
    high_var_genes = gene_data.columns[gene_data.var() > 0]
    selected_gene_data = gene_data[high_var_genes]
    target_median = selected_gene_data.median(axis=1).median(axis=0)
    scale_factors = target_median / selected_gene_data.median(axis=1)
    expression_data = selected_gene_data.multiply(scale_factors, axis=0)
    print(f"[Filtered Genes] {expression_data.shape[1]}")

    # Perform DEG analysis
    print(f"[Differential Expression] ...", end='\r')
    group1_data = expression_data[metadata[group_column].isin(group1)]
    group2_data = expression_data[metadata[group_column].isin(group2)]
    
    if len(group1_data) == 0 or len(group2_data) == 0:
        print(f"Error: One or both groups have no samples. Skipping differential expression analysis.")
        return None
    
    pvalues = []
    log2_fold_changes = []
    genes = []
    
    for gene in expression_data.columns:
        if group1_data[gene].nunique() == 1 and group2_data[gene].nunique() == 1 and group1_data[gene].iloc[0] == group2_data[gene].iloc[0]:
            continue
        
        t_stat, p_value = stats.ttest_ind(group1_data[gene], group2_data[gene])
        if np.isnan(p_value):
            continue
        pvalues.append(p_value)
        genes.append(gene)
        
        mean1 = group1_data[gene].mean()
        mean2 = group2_data[gene].mean()
        if mean1 == 0 and mean2 == 0:
            log2_fold_changes.append(0)
        elif mean1 == 0:
            log2_fold_changes.append(np.inf)
        elif mean2 == 0:
            log2_fold_changes.append(-np.inf)
        else:
            log2_fold_changes.append(np.log2(mean2 / mean1))
    
    _, adjusted_pvalues, _, _ = multipletests(pvalues, method='fdr_bh')
    
    results_df = pd.DataFrame({
        'gene': genes,
        'log2_fold_change': log2_fold_changes,
        'p_value': pvalues,
        'adjusted_pvalue': adjusted_pvalues,
    })

    # Count genes
    log2_fc_threshold = 1
    p_value_threshold = 0.05
    down_regulated = sum((results_df['log2_fold_change'] < -log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold))
    up_regulated = sum((results_df['log2_fold_change'] > log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold))

    # Save results
    results_df = results_df.sort_values('adjusted_pvalue')
    output_dir = os.path.dirname(summarize_file)
    results_file = os.path.join(output_dir, f'DEG_{"+".join(group1)}_vs_{"+".join(group2)}.csv')
    results_df.to_csv(results_file, index=False)
    print(f"[Differential Expression] Down: {down_regulated} / Up: {up_regulated}")

    # Generate volcano plot
    print(f"[Volcano Plot] ...", end='\r')
    plt.figure(figsize=(10, 10))
    plt.style.use('ggplot')

    results_df['color'] = 'grey'
    results_df.loc[(results_df['log2_fold_change'] > log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold), 'color'] = 'red'
    results_df.loc[(results_df['log2_fold_change'] < -log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold), 'color'] = 'blue'
    
    plt.scatter(results_df['log2_fold_change'], -np.log10(results_df['adjusted_pvalue']), c=results_df['color'], alpha=0.5)
    
    plt.axvline(x=log2_fc_threshold, color='gray', linestyle='--')
    plt.axvline(x=-log2_fc_threshold, color='gray', linestyle='--')
    plt.axhline(y=-np.log10(p_value_threshold), color='gray', linestyle='--')
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-value')
    
    volcano_file = os.path.join(output_dir, f'volcano_{"+".join(group1)}_vs_{"+".join(group2)}.pdf')
    plt.savefig(volcano_file, format='pdf', dpi=600, bbox_inches='tight')
    volcano_file = os.path.join(output_dir, f'volcano_{"+".join(group1)}_vs_{"+".join(group2)}.png')
    plt.savefig(volcano_file, format='png', dpi=600, bbox_inches='tight')
    plt.close()
    print("[Volcano Plot] Complete")

    # Generate heatmap
    print(f"[Generate Heatmap] ...", end='\r')
    filtered_metadata = metadata[metadata[group_column].isin(group1 + group2)]
    filtered_expression_data = expression_data.loc[filtered_metadata.index]

    down_genes = results_df[(results_df['adjusted_pvalue'] < p_value_threshold) & (results_df['log2_fold_change'] < -log2_fc_threshold)]['gene']
    up_genes = results_df[(results_df['adjusted_pvalue'] < p_value_threshold) & (results_df['log2_fold_change'] > log2_fc_threshold)]['gene']
    
    significant_genes = pd.concat([down_genes, up_genes])
    selected_gene_data = filtered_expression_data[significant_genes]
    
    # Cluster genes
    gene_order = []
    for gene_group in [down_genes, up_genes]:
        gene_data = selected_gene_data[gene_group]
        if len(gene_data.columns) > 1:
            gene_dist = pdist(gene_data.T, metric='correlation')
            gene_linkage = hierarchy.linkage(gene_dist, method='ward')
            gene_order.extend(gene_data.columns[hierarchy.leaves_list(gene_linkage)])
        else:
            gene_order.extend(gene_data.columns)

    # Cluster samples
    group_order = []
    for group in [group1, group2]:
        group_data = selected_gene_data[filtered_metadata[group_column].isin(group)]
        if len(group_data) > 1:
            group_dist = pdist(group_data, metric='correlation')
            group_linkage = hierarchy.linkage(group_dist, method='ward')
            group_order.extend(group_data.index[hierarchy.leaves_list(group_linkage)])
        else:
            group_order.extend(group_data.index)

    # Reorder the data
    data_ordered = selected_gene_data.loc[group_order, gene_order]
    
    # Calculate z-scores for visualization
    scaler = StandardScaler()
    z_scores = pd.DataFrame(scaler.fit_transform(data_ordered), 
                            columns=data_ordered.columns, 
                            index=data_ordered.index)

    fig = plt.figure(figsize=(11.3, 10.3))
    gs = fig.add_gridspec(2, 3, width_ratios=[0.3, 10, 1], height_ratios=[0.3, 10],
                          left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.02, hspace=0.02)
    
    # Heatmap
    ax_heatmap = fig.add_subplot(gs[1, 1])
    vmin = np.percentile(z_scores.values, 1)
    vmax = np.percentile(z_scores.values, 99)
    sns.heatmap(z_scores.T, cmap='seismic', center=0, vmin=vmin, vmax=vmax,
                xticklabels=False, yticklabels=False, cbar=False, ax=ax_heatmap)
    
    # Group indicator (top)
    ax_groups = fig.add_subplot(gs[0, 1], sharex=ax_heatmap)
    ax_groups.set_ylim(0, 1)
    
    cmap = plt.get_cmap('seismic')
    color_map = {0: cmap(0.9), 1: cmap(0.1)}
    start = 0
    for i, group in enumerate([group1, group2]):
        group_samples = filtered_metadata[filtered_metadata[group_column].isin(group)]
        width = len(group_samples)
        ax_groups.add_patch(Rectangle((start, 0), width, 1, facecolor=color_map[i], edgecolor='none'))
        ax_groups.text(start + width/2, 0.5, "+".join(group), ha='center', va='center', color='white', fontweight='semibold')
        start += width
    
    ax_groups.set_xlim(0, len(data_ordered))
    ax_groups.axis('off')
    
    # Gene regulation indicator (left)
    ax_regulation = fig.add_subplot(gs[1, 0], sharey=ax_heatmap)
    ax_regulation.set_xlim(0, 1)
    down_height = len(down_genes)
    up_height = len(up_genes)
    ax_regulation.add_patch(Rectangle((0, 0), 1, down_height, facecolor=color_map[0], edgecolor='none'))
    ax_regulation.add_patch(Rectangle((0, down_height), 1, up_height, facecolor=color_map[1], edgecolor='none'))
    ax_regulation.text(0.5, down_height/2, 'Down', ha='center', va='center', color='white', fontweight='semibold', rotation=90)
    ax_regulation.text(0.5, down_height + up_height/2, 'Up', ha='center', va='center', color='white', fontweight='semibold', rotation=90)
    ax_regulation.set_ylim(down_height + up_height, 0)
    ax_regulation.axis('off')
    
    # Legend and colorbar
    ax_legend = fig.add_subplot(gs[1, 2])

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
    
    legend = ax_legend.legend(handles=cbar_elements, loc='center', 
                              ncol=1, handlelength=1, handleheight=1, 
                              handletextpad=0.5, columnspacing=0.5, labelspacing=0.0)
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_facecolor('none')
    ax_legend.axis('off')
    
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Save the plot
    heatmap_file = os.path.join(output_dir, f'heatmap_{"+".join(group1)}_vs_{"+".join(group2)}.pdf')
    plt.savefig(heatmap_file, format='pdf', dpi=600, bbox_inches='tight')
    heatmap_file = os.path.join(output_dir, f'heatmap_{"+".join(group1)}_vs_{"+".join(group2)}.png')
    plt.savefig(heatmap_file, format='png', dpi=600, bbox_inches='tight')
    plt.close()
    print("[Generate Heatmap] Complete")


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python3 DEG.py <cohort/summarize.csv> <group_column> <group1a> <group1b> ... -- <group2a> <group2b> ...")
        sys.exit(1)
    
    summarize_file = sys.argv[1]
    group_column = sys.argv[2]
    
    try:
        separator_index = sys.argv.index('--')
        group1 = sys.argv[3:separator_index]
        group2 = sys.argv[separator_index+1:]
    except ValueError:
        print("Error: Missing '--' separator between groups")
        sys.exit(1)
    
    if not group1 or not group2:
        print("Error: Both groups must contain at least one class")
        sys.exit(1)
    
    main(summarize_file, group_column, group1, group2)