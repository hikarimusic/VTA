# -------------------------

gene_threshold = 1
gene_normalize = True

deg_multiple_test_correction = 'fdr_bh'
deg_log2_fc_threshold = 1
deg_p_value_threshold = 0.05

count_plot_generate = False
count_plot_format = 'png'
count_plot_size = (3.5, 3.5)
count_plot_fontsize = 6
count_plot_dpi = 600

volcano_plot_format = 'png'
volcano_plot_size = (3.5, 3.5)
volcano_plot_fontsize = 6
volcano_plot_dotsize = 6
volcano_plot_dotalpha = 0.5
volcano_plot_up_color = 'red'
volcano_plot_down_color = 'blue'
volcano_plot_other_color = 'grey'
volcano_plot_line = True
volcano_plot_dpi = 600

strip_plot_format = 'png'
strip_plot_size = (3.5, 7.0)
strip_plot_fontsize = 6
strip_plot_dotsize = 3
strip_plot_dotalpha = 0.5
strip_plot_color = 'seismic'
strip_plot_number = 50
strip_plot_dpi = 600

heatmap_format = 'png'
heatmap_size = (3.5, 3.5)
heatmap_color = 'seismic'
heatmap_fontsize = 6
heatmap_gene_metric = 'correlation'
heatmap_gene_method = 'ward'
heatmap_dpi = 600

# -------------------------

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

def generate_count_plot(metadata, group_column, group1, group2, output_dir):
    print("[Count Plot] ...                 ", end='\r')

    # Combined group
    filtered_metadata = metadata[metadata[group_column].isin(group1 + group2)].copy()
    filtered_metadata['Group'] = filtered_metadata[group_column].apply(
        lambda x: '+'.join(group1) if x in group1 else '+'.join(group2)
    )
    
    # Iterate group columns
    for col in metadata.columns:
        if col != group_column and col != 'Group':
            unique_values = filtered_metadata[col].nunique()
            if unique_values >=2 and unique_values <= 10:
                # Chi-square test
                contingency_table = pd.crosstab(filtered_metadata['Group'], filtered_metadata[col])
                p_value = "None"
                if (contingency_table > 0).all().all():
                    chi2, p_value = stats.chi2_contingency(contingency_table)[:2]
                    if p_value < 0.001:
                        p_value = f'{p_value:.2e}'
                    else:
                        p_value = f'{p_value:.3f}'

                # Count plot
                plt.style.use('ggplot')
                plt.figure(figsize=count_plot_size)
                sns.countplot(data=filtered_metadata, x='Group', hue=col)

                plt.xticks(fontsize=count_plot_fontsize)
                plt.yticks(fontsize=count_plot_fontsize)
                plt.xlabel('Group', fontsize=count_plot_fontsize)
                plt.ylabel('Count', fontsize=count_plot_fontsize)
                plt.legend(fontsize=count_plot_fontsize)

                count_plot_file = os.path.join(output_dir, f'DEA_count_{"+".join(group1)}_vs_{"+".join(group2)}_{col}.' + count_plot_format)
                plt.savefig(count_plot_file, format=count_plot_format, dpi=count_plot_dpi, bbox_inches='tight')
                plt.close()

                print(f"[Chi-Square] Group / {col}: {p_value}                 ")
    
    print("[Count Plot] Complete                 ")

def generate_volcano_plot(results_df, group1, group2, output_dir, log2_fc_threshold, p_value_threshold):
    print("[Volcano Plot] ...                 ", end='\r')

    # Up and down genes
    results_df['color'] = volcano_plot_other_color
    results_df.loc[(results_df['log2_fold_change'] > log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold), 'color'] = volcano_plot_up_color
    results_df.loc[(results_df['log2_fold_change'] < -log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold), 'color'] = volcano_plot_down_color

    # Volcano plot
    plt.style.use('ggplot')
    plt.figure(figsize=volcano_plot_size)
    plt.scatter(results_df['log2_fold_change'], -np.log10(results_df['adjusted_pvalue']), c=results_df['color'], alpha=volcano_plot_dotalpha, s=volcano_plot_dotsize)
    
    if volcano_plot_line == True:
        plt.axvline(x=log2_fc_threshold, color='gray', linestyle='--')
        plt.axvline(x=-log2_fc_threshold, color='gray', linestyle='--')
        plt.axhline(y=-np.log10(p_value_threshold), color='gray', linestyle='--')
    
    plt.xticks(fontsize=volcano_plot_fontsize)
    plt.yticks(fontsize=volcano_plot_fontsize)
    plt.xlabel('Log2 Fold Change', fontsize=volcano_plot_fontsize)
    plt.ylabel('-Log10 Adjusted P-value', fontsize=volcano_plot_fontsize)
    
    volcano_file = os.path.join(output_dir, f'DEA_volcano_{"+".join(group1)}_vs_{"+".join(group2)}.' + volcano_plot_format)
    plt.savefig(volcano_file, format=volcano_plot_format, dpi=volcano_plot_dpi, bbox_inches='tight')
    plt.close()

    print("[Volcano Plot] Complete                 ")

def generate_strip_plot(results_df, expression_data, metadata, group_column, group1, group2, output_dir, log2_fc_threshold, p_value_threshold):
    print("[Strip Plot] ...                 ", end='\r')
    
    # Pre-process
    sig_genes = results_df[results_df['adjusted_pvalue'] < p_value_threshold].copy()
    up_genes = sig_genes[sig_genes['log2_fold_change'] > log2_fc_threshold].nsmallest(strip_plot_number, 'p_value')['gene']
    down_genes = sig_genes[sig_genes['log2_fold_change'] < -log2_fc_threshold].nsmallest(strip_plot_number, 'p_value')['gene']
    
    filtered_metadata = metadata[metadata[group_column].isin(group1 + group2)].copy()
    filtered_metadata['Group'] = filtered_metadata[group_column].apply(
        lambda x: '+'.join(group1) if x in group1 else '+'.join(group2)
    )
    
    scaler = StandardScaler()
    z_scores = pd.DataFrame(
        scaler.fit_transform(expression_data),
        columns=expression_data.columns,
        index=expression_data.index
    )
    
    # Strip plot
    cmap = plt.get_cmap(strip_plot_color)
    color_map = {0: cmap(0.1), 1: cmap(0.9)}
    palette = {'+'.join(group1): color_map[0], '+'.join(group2): color_map[1]}
    for gene_set, title_prefix in [(up_genes, 'up'), (down_genes, 'down')]:
        if len(gene_set) == 0:
            continue
        
        plot_data = pd.DataFrame()
        for gene in gene_set:
            gene_data = z_scores[gene]
            temp_df = pd.DataFrame({
                'Gene': gene,
                'Z-score': gene_data,
                'Group': filtered_metadata['Group']
            })
            plot_data = pd.concat([plot_data, temp_df])

        plt.style.use('ggplot')
        plt.figure(figsize=strip_plot_size)
        sns.stripplot(data=plot_data, x='Z-score', y='Gene', hue='Group', palette=palette, alpha=strip_plot_dotalpha, s=strip_plot_dotsize)
        
        plt.xticks(fontsize=strip_plot_fontsize)
        plt.yticks(fontsize=strip_plot_fontsize)
        plt.xlabel('Expression Z-score', fontsize=strip_plot_fontsize)
        plt.ylabel('Gene', fontsize=strip_plot_fontsize)
        plt.legend(loc='upper right', fontsize=strip_plot_fontsize)
        
        strip_plot_file = os.path.join(output_dir, f'DEA_strip_{"+".join(group1)}_vs_{"+".join(group2)}_{title_prefix}.' + strip_plot_format)
        plt.savefig(strip_plot_file, format=strip_plot_format, dpi=strip_plot_dpi, bbox_inches='tight')
        plt.close()
    
    print("[Strip Plot] Complete                 ")

def generate_heatmap(results_df, expression_data, metadata, group_column, group1, group2, output_dir, log2_fc_threshold, p_value_threshold):
    print("[Heatmap] ...                 ", end='\r')

    # Pre-process
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
            gene_dist = pdist(gene_data.T, metric=heatmap_gene_metric)
            gene_linkage = hierarchy.linkage(gene_dist, method=heatmap_gene_method)
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

    data_ordered = selected_gene_data.loc[group_order, gene_order]
    
    # Calculate z-scores
    scaler = StandardScaler()
    z_scores = pd.DataFrame(scaler.fit_transform(data_ordered), 
                            columns=data_ordered.columns, 
                            index=data_ordered.index)

    # Start plotting
    plt.style.use('seaborn-v0_8-whitegrid')
    fig = plt.figure(figsize=heatmap_size)
    gs = fig.add_gridspec(2, 3, width_ratios=[0.4, 10, 1], height_ratios=[0.4, 10],
                          left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.02, hspace=0.02)
    
    # Heatmap
    ax_heatmap = fig.add_subplot(gs[1, 1])
    vmin = np.percentile(z_scores.values, 1)
    vmax = np.percentile(z_scores.values, 99)
    sns.heatmap(z_scores.T, cmap=heatmap_color, center=0, vmin=vmin, vmax=vmax,
                xticklabels=False, yticklabels=False, cbar=False, ax=ax_heatmap)
    
    # Group indicator (top)
    ax_groups = fig.add_subplot(gs[0, 1], sharex=ax_heatmap)
    ax_groups.set_ylim(0, 1)
    
    cmap = plt.get_cmap(heatmap_color)
    color_map = {0: cmap(0.1), 1: cmap(0.9)}
    start = 0
    for i, group in enumerate([group1, group2]):
        group_samples = filtered_metadata[filtered_metadata[group_column].isin(group)]
        width = len(group_samples)
        ax_groups.add_patch(Rectangle((start, 0), width, 1, facecolor=color_map[i], edgecolor='none'))
        ax_groups.text(start + width/2, 0.5, "+".join(group), ha='center', va='center', color='white', fontsize=heatmap_fontsize, fontweight='semibold')
        start += width
    ax_groups.set_xlim(0, len(data_ordered))
    ax_groups.axis('off')
    
    # Regulation indicator (left)
    ax_regulation = fig.add_subplot(gs[1, 0], sharey=ax_heatmap)
    ax_regulation.set_xlim(0, 1)
    down_height = len(down_genes)
    up_height = len(up_genes)
    ax_regulation.add_patch(Rectangle((0, 0), 1, down_height, facecolor=color_map[0], edgecolor='none'))
    ax_regulation.add_patch(Rectangle((0, down_height), 1, up_height, facecolor=color_map[1], edgecolor='none'))
    ax_regulation.text(0.5, down_height/2, 'Down', ha='center', va='center', color='white', fontsize=heatmap_fontsize, fontweight='semibold', rotation=90)
    ax_regulation.text(0.5, down_height + up_height/2, 'Up', ha='center', va='center', color='white', fontsize=heatmap_fontsize, fontweight='semibold', rotation=90)
    ax_regulation.set_ylim(down_height + up_height, 0)
    ax_regulation.axis('off')
    
    # Legend and colorbar
    ax_legend = fig.add_subplot(gs[1, 2])
    cmap = plt.get_cmap(heatmap_color)
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
                              handletextpad=0.5, columnspacing=0.5, labelspacing=0.0,
                              prop={'size': heatmap_fontsize})
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_facecolor('none')
    ax_legend.axis('off')
    
    # Save the plot
    heatmap_file = os.path.join(output_dir, f'DEA_heatmap_{"+".join(group1)}_vs_{"+".join(group2)}.' + heatmap_format)
    plt.savefig(heatmap_file, format=heatmap_format, dpi=heatmap_dpi, bbox_inches='tight')
    plt.close()

    print("[Heatmap] Complete                 ")

def DEG(summarize_file, group_column, group1, group2):
    # Read data
    print(f"[Read Data] ...", end='\r')
    df = pd.read_csv(summarize_file, low_memory=False)
    df = df.fillna('~')
    start_gene_index = df.columns.get_loc('START_GENE')
    metadata = df.iloc[:, :start_gene_index + 1]
    gene_data = df.iloc[:, start_gene_index + 1:]
    print("[Read Data] Complete                 ")

    # Filter and normalize gene
    print(f"[Filter Genes] ...", end='\r')
    expressed_genes = gene_data.columns[gene_data.mean() > gene_threshold]
    gene_data = gene_data[expressed_genes]
    high_var_genes = gene_data.columns[gene_data.var() > 0]
    selected_gene_data = gene_data[high_var_genes]
    target_median = selected_gene_data.median(axis=1).median(axis=0)
    scale_factors = target_median / selected_gene_data.median(axis=1)
    if gene_normalize == False:
        scale_factors = 1
    expression_data = selected_gene_data.multiply(scale_factors, axis=0)
    print(f"[Filter Genes] {expression_data.shape[1]}                 ")
     
    # Perform DEG analysis
    print(f"[Find DEG] ...", end='\r')
    group1_data = expression_data[metadata[group_column].isin(group1)]
    group2_data = expression_data[metadata[group_column].isin(group2)]
    
    genes = []
    pvalues = []
    log2_fold_changes = []
    for gene in expression_data.columns:
        mean1 = group1_data[gene].mean()
        mean2 = group2_data[gene].mean()
        if mean1 < 1e-9 or mean2 < 1e-9:
            continue
        log2_fold_change = np.log2(mean2 / mean1)
        _, p_value = stats.ttest_ind(group1_data[gene], group2_data[gene])
        genes.append(gene)
        pvalues.append(p_value)
        log2_fold_changes.append(log2_fold_change)
    
    # Result df
    _, adjusted_pvalues, _, _ = multipletests(pvalues, method=deg_multiple_test_correction)
    results_df = pd.DataFrame({
        'gene': genes,
        'log2_fold_change': log2_fold_changes,
        'p_value': pvalues,
        'adjusted_pvalue': adjusted_pvalues,
    })

    # Count significant genes
    log2_fc_threshold = deg_log2_fc_threshold
    p_value_threshold = deg_p_value_threshold
    down_regulated = sum((results_df['log2_fold_change'] < -log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold))
    up_regulated = sum((results_df['log2_fold_change'] > log2_fc_threshold) & (results_df['adjusted_pvalue'] < p_value_threshold))    
    print(f"[Find DEG] Down: {down_regulated} / Up: {up_regulated}                 ")

    # Save results
    print("[Save Genes] ...", end='\r')
    results_df = results_df.sort_values('p_value')

    output_dir = os.path.dirname(summarize_file)
    results_file = os.path.join(output_dir, f'DEA_genes_{"+".join(group1)}_vs_{"+".join(group2)}.csv')
    results_df.to_csv(results_file, index=False)
    print(f"[Save Genes] Complete                 ")

    print("[Create Plots] ...", end='\r')

    # Count plot
    if count_plot_generate == True:
        generate_count_plot(metadata, group_column, group1, group2, output_dir)
    
    # Volcano plot
    generate_volcano_plot(results_df, group1, group2, output_dir, log2_fc_threshold, p_value_threshold)

    # Strip plots
    generate_strip_plot(results_df, expression_data, metadata, group_column, group1, group2, output_dir, log2_fc_threshold, p_value_threshold)

    # Heatmap
    generate_heatmap(results_df, expression_data, metadata, group_column, group1, group2, output_dir, log2_fc_threshold, p_value_threshold)

    print(f"[Create Plots] Complete                 ")


if __name__ == "__main__":
    '''Command: python3 DEA.py <cohort/summary.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... '''
    summarize_file = sys.argv[1]
    group_column = sys.argv[2]
    separator_index = sys.argv.index('--')
    group1 = sys.argv[3:separator_index]
    group2 = sys.argv[separator_index+1:]
    
    DEG(summarize_file, group_column, group1, group2)