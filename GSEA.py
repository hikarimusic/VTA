# -------------------------

gene_threshold = 1
gene_normalize = True

ks_multiple_test_correction = 'fdr_bh'
ks_position_threshold = (0.40, 0.60)
ks_p_value_threshold = 0.001

gsea_plot_format = 'png'
gsea_plot_size = (3.5, 3.5)
gsea_plot_dpi = 600
gsea_plot_fontsize = 6
gsea_plot_linewidth = 1
gsea_plot_linealpha = 0.7
gsea_plot_hit_prop = 0.05
gsea_plot_hitwidth = 1
gsea_plot_hitalpha = 0.2
gsea_plot_up_color = 'red'
gsea_plot_down_color = 'blue'
gsea_plot_other_color = 'gray'
gsea_plot_hit_color = 'black'

bar_plot_format = 'png'
bar_plot_size = (3.5, 7)
bar_plot_dpi = 600
bar_plot_fontsize = 6
bar_plot_alpha = 0.7
bar_plot_up_color = 'red'
bar_plot_down_color = 'blue'

# -------------------------

import sys
import os
import shutil
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_correlation(group1_data, group2_data):
    if group1_data.empty or group2_data.empty:
        return 0

    mean1 = group1_data.mean()
    mean2 = group2_data.mean()
    log2_fold_change = np.log2((mean2 + 1e-9) / (mean1 + 1e-9))
    _, p_value = stats.ttest_ind(group1_data, group2_data)
    correlation = np.sign(log2_fold_change) * -np.log10(p_value)

    return correlation

def calculate_ks_test(ranked_genes, gene_set):
    # Hit genes
    hit_indices = [i for i, gene in enumerate(ranked_genes) if gene in gene_set]
    if not hit_indices:
        return 0, 1, 0
    
    # Kolmogorov-Smirnov test

    results = stats.kstest(hit_indices, lambda x: stats.uniform.cdf(x, loc=0, scale=len(ranked_genes)-1))
    
    enrichment = results.statistic * results.statistic_sign
    position = results.statistic_location / len(ranked_genes)
    p_value = results.pvalue

    return enrichment, position, p_value

def calculate_enrichment_score(ranked_genes, gene_set):
    N = len(ranked_genes)
    hits = [1 if gene in gene_set else 0 for gene in ranked_genes]
    hit_indices = np.where(np.array(hits) == 1)[0]
    
    # Running enrichment score
    running_sum = np.zeros(N)
    total_hits = sum(hits)

    hit_weight = 1.0 / total_hits if total_hits > 0 else 0
    miss_weight = 1.0 / N
    curr_sum = 0
    for i in range(N):
        if hits[i] == 1:
            curr_sum += hit_weight
        curr_sum -= miss_weight
        running_sum[i] = curr_sum
    
    return running_sum, hit_indices

def plot_gsea_plot(ranked_genes, gene_set, gene_set_name, p_adjust, correlations, plots_dir, direction='up'):
    running_sum, hit_indices = calculate_enrichment_score(ranked_genes, gene_set)
    
    # Find key x-coordinates
    max_es_idx = np.argmax(running_sum) if direction == 'up' else np.argmin(running_sum)
    correlation_values = [correlations[gene] for gene in ranked_genes]

    # Create GSEA plot
    plt.style.use('ggplot')
    plt.figure(figsize=gsea_plot_size)
    
    x = np.arange(len(ranked_genes))
    plt.plot(x, running_sum, color=gsea_plot_up_color if direction == 'up' else gsea_plot_down_color, linewidth=gsea_plot_linewidth)
    plt.axhline(y=0, color=gsea_plot_other_color, linestyle='--', alpha=gsea_plot_linealpha, linewidth=gsea_plot_linewidth)
    plt.axvline(x=max_es_idx, color=gsea_plot_up_color if direction == 'up' else gsea_plot_down_color, linestyle='--', alpha=gsea_plot_linealpha, linewidth=gsea_plot_linewidth)

    # Mark hits
    v_max = np.max(running_sum)
    v_min = np.min(running_sum)
    v_hit = (v_max - v_min) * (1 if np.abs(v_max) > np.abs(v_min) else -1) * gsea_plot_hit_prop
    plt.vlines(hit_indices, 0, v_hit, color=gsea_plot_hit_color, alpha=gsea_plot_hitalpha, linewidth=gsea_plot_hitwidth)
    
    # Customize
    plt.xticks(fontsize=gsea_plot_fontsize)
    plt.yticks(fontsize=gsea_plot_fontsize)
    plt.ylabel('Enrichment Score', fontsize=gsea_plot_fontsize)
    plt.xlabel('Rank of Genes', fontsize=gsea_plot_fontsize)
    plt.xlim(0, len(ranked_genes))
    plt.ylim(np.min(running_sum) - 0.02, np.max(running_sum) + 0.02)
    
    # Text box
    gene_set_text = '_'.join(gene_set_name.split('_')[1:])
    if p_adjust < 0.001:
        p_value_text = f'p-adjust = {p_adjust:.2e}'
    else:
        p_value_text = f'p-adjust = {p_adjust:.3f}'
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, fc='none', ec='none', label=gene_set_text),
        plt.Rectangle((0, 0), 1, 1, fc='none', ec='none', label=p_value_text)
    ]
    legend = plt.legend(handles=legend_elements, 
                        loc='upper right' if direction == 'up' else 'lower left',
                        frameon=False,
                        fontsize=gsea_plot_fontsize,
                        handlelength=0)

    # Save GSEA plot
    plt.tight_layout()
    gsea_output = os.path.join(plots_dir, f'{direction}_{gene_set_name}.' + gsea_plot_format)
    plt.savefig(gsea_output, dpi=gsea_plot_dpi, bbox_inches='tight')
    plt.close()


def plot_bar_plot(ranked_genes, gene_set, gene_set_name, p_adjust, correlations, plots_dir, direction='up'):
    running_sum, hit_indices = calculate_enrichment_score(ranked_genes, gene_set)
    max_es_idx = np.argmax(running_sum) if direction == 'up' else np.argmin(running_sum)
    
    # Create bar plot
    hit_genes = [ranked_genes[i] for i in hit_indices]
    if direction == 'up':
        selected_hits = [gene for i, gene in enumerate(hit_genes) if hit_indices[i] < max_es_idx]
    else:
        selected_hits = [gene for i, gene in enumerate(hit_genes) if hit_indices[i] > max_es_idx]
        selected_hits.reverse()
   
    if selected_hits:
        p_values = [abs(correlations[gene]) for gene in selected_hits]
        plt.figure(figsize=bar_plot_size)
        bars = plt.barh(range(len(selected_hits)), p_values, 
                        color=bar_plot_up_color if direction == 'up' else bar_plot_down_color, alpha=bar_plot_alpha)
       
        # Customize
        plt.yticks(range(len(selected_hits)), selected_hits, rotation=0, fontsize=bar_plot_fontsize)
        plt.xticks(fontsize=bar_plot_fontsize)
        plt.xlabel('-Log10 P-value', fontsize=bar_plot_fontsize)
        plt.ylabel('Leading Genes', fontsize=bar_plot_fontsize)
        plt.ylim(len(selected_hits), -1)
        
        # Text box
        gene_set_text = '_'.join(gene_set_name.split('_')[1:])
        gene_cnt_text = f'Gene Count = {len(selected_hits)}'
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, fc='none', ec='none', label=gene_set_text),
            plt.Rectangle((0, 0), 1, 1, fc='none', ec='none', label=gene_cnt_text)
        ]
        legend = plt.legend(handles=legend_elements,
                            loc='lower right',
                            frameon=False,
                            fontsize=bar_plot_fontsize,
                            handlelength=0)
       
        # Save bar plot
        plt.tight_layout()
        bar_output = os.path.join(plots_dir, f'{direction}_{gene_set_name}.' + bar_plot_format)
        plt.savefig(bar_output, dpi=bar_plot_dpi, bbox_inches='tight')
        plt.close()
    

def GSEA(summarize_file, group_column, group1, group2, gmt_file):    
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

    # Load GeneSets
    print(f"[Load GeneSets] ...", end='\r')
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            gene_set_name = parts[0]
            genes = set(parts[2:])
            gene_sets[gene_set_name] = genes
    print(f"[Load GeneSets] {len(gene_sets)}                 ")   

    # Perform GSEA with KS test
    print("[KS Test] ...", end='\r')
    phenotypes = metadata[group_column].map(lambda x: 1 if x in group1 else (-1 if x in group2 else 0))
    expression_data = expression_data.loc[phenotypes != 0]
    phenotypes = phenotypes[phenotypes != 0]
    genes = expression_data.columns.tolist()

    group1_data = expression_data[phenotypes == 1]
    group2_data = expression_data[phenotypes == -1]
    correlations = {gene: calculate_correlation(group1_data[gene], group2_data[gene]) for gene in genes}
    ranked_genes = sorted(correlations, key=correlations.get, reverse=True)

    results = []
    for i, (gene_set_name, gene_set) in enumerate(gene_sets.items(), 1):
        enrichment, position, p_value = calculate_ks_test(ranked_genes, gene_set)
        results.append((gene_set_name, enrichment, position, p_value))
    
    # Adjust p-values
    _, adjusted_pvalues, _, _ = multipletests([r[3] for r in results], method=ks_multiple_test_correction)
    results = [(r[0], r[1], r[2], r[3], adj_p) for r, adj_p in zip(results, adjusted_pvalues)]
    results_df = pd.DataFrame(results, columns=['gene_set', 'enrichment_score', 'position', 'p_value', 'adjusted_pvalue'])
    up_regulated = sum((results_df['enrichment_score'] > 0) & 
                       (results_df['position'] < ks_position_threshold[0]) & 
                       (results_df['adjusted_pvalue'] < ks_p_value_threshold))  
    down_regulated = sum((results_df['enrichment_score'] < 0) & 
                         (results_df['position'] > ks_position_threshold[1]) & 
                         (results_df['adjusted_pvalue'] < ks_p_value_threshold))  
    print(f"[KS Test] Up: {up_regulated} / Down: {down_regulated}                 ")

    # Save results
    print("[Save GeneSets] ...", end='\r')
    geneset_name = os.path.splitext(os.path.basename(gmt_file))[0]
    output_dir = os.path.dirname(summarize_file)

    results_df = results_df.sort_values('p_value')
    results_file = os.path.join(output_dir, f'GSEA_genesets_{"+".join(group1)}_vs_{"+".join(group2)}_{geneset_name}.csv')
    results_df.to_csv(results_file, index=False)
    print("[Save GeneSets] Complete                 ")

    # GSEA plots
    print("[GSEA Plots] ...", end='\r')
    plots_dir = os.path.join(output_dir, f'GSEA_gsea_{"+".join(group1)}_vs_{"+".join(group2)}_{geneset_name}')
    if os.path.exists(plots_dir):
        shutil.rmtree(plots_dir)
    os.makedirs(plots_dir)

    geneset_cnt = 0
    for _, row in results_df.iterrows():
        gene_set_name = row['gene_set']
        enrichment_score = row['enrichment_score']
        position = row['position']
        p_adjust = row['adjusted_pvalue']
        if p_adjust < ks_p_value_threshold:
            if (enrichment_score > 0) and (position < ks_position_threshold[0]):
                plot_gsea_plot(ranked_genes, gene_sets[gene_set_name], gene_set_name, p_adjust, 
                                               correlations, plots_dir, direction='up')
                geneset_cnt += 1
            elif (enrichment_score < 0) and (position > ks_position_threshold[1]):
                plot_gsea_plot(ranked_genes, gene_sets[gene_set_name], gene_set_name, p_adjust, 
                                               correlations, plots_dir, direction='down')
                geneset_cnt += 1
            print(f"[GSEA Plots] {geneset_cnt} / {up_regulated+down_regulated}", end='\r')

    print("[GSEA Plots] Complete                 ")
    
    # Bar plots
    print("[Bar Plots] ...", end='\r')
    plots_dir = os.path.join(output_dir, f'GSEA_bar_{"+".join(group1)}_vs_{"+".join(group2)}_{geneset_name}')
    if os.path.exists(plots_dir):
        shutil.rmtree(plots_dir)
    os.makedirs(plots_dir)

    geneset_cnt = 0
    for _, row in results_df.iterrows():
        gene_set_name = row['gene_set']
        enrichment_score = row['enrichment_score']
        position = row['position']
        p_adjust = row['adjusted_pvalue']
        if p_adjust < ks_p_value_threshold:
            if (enrichment_score > 0) and (position < ks_position_threshold[0]):
                plot_bar_plot(ranked_genes, gene_sets[gene_set_name], gene_set_name, p_adjust, 
                                               correlations, plots_dir, direction='up')
                geneset_cnt += 1
            elif (enrichment_score < 0) and (position > ks_position_threshold[1]):
                plot_bar_plot(ranked_genes, gene_sets[gene_set_name], gene_set_name, p_adjust, 
                                               correlations, plots_dir, direction='down')
                geneset_cnt += 1
            print(f"[Bar Plots] {geneset_cnt} / {up_regulated+down_regulated}", end='\r')

    print("[Bar Plots] Complete                 ")
    

if __name__ == "__main__":
    '''Command: python3 GSEA.py <cohort/summary.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... <geneset.gmt>'''
    summarize_file = sys.argv[1]
    group_column = sys.argv[2]
    separator_index = sys.argv.index('--')
    group1 = sys.argv[3:separator_index]
    group2 = sys.argv[separator_index+1:-1]
    gmt_file = sys.argv[-1]
    
    GSEA(summarize_file, group_column, group1, group2, gmt_file)