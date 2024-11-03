import sys
import os
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
    log2_fold_change = np.log2((mean2 + 1e-8) / (mean1 + 1e-8))
    _, p_value = stats.ttest_ind(group1_data, group2_data)
    correlation = np.sign(log2_fold_change) * -np.log10(p_value)

    return correlation

def calculate_ks_test(ranked_genes, gene_set, direction='up'):
    # Hit genes
    hit_indices = [i for i, gene in enumerate(ranked_genes) if gene in gene_set]
    if not hit_indices:
        return 0, 1, 0
    
    # Kolmogorov-Smirnov test
    ks_stat, p_value = stats.kstest(hit_indices, 
                                    lambda x: stats.uniform.cdf(x, loc=0, scale=len(ranked_genes)-1), 
                                    alternative='greater' if direction == 'up' else 'less')
    return ks_stat, p_value, 1 if direction == 'up' else -1

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

def plot_gsea_bar(ranked_genes, gene_set, gene_set_name, p_adjust, correlations, plots_dir, direction='up'):
    running_sum, hit_indices = calculate_enrichment_score(ranked_genes, gene_set)
    
    # Find key x-coordinates
    max_es_idx = np.argmax(running_sum) if direction == 'up' else np.argmin(running_sum)
    correlation_values = [correlations[gene] for gene in ranked_genes]
    zero_cross_indices = np.where(np.diff(np.signbit(correlation_values)))[0]
    zero_cross_idx = zero_cross_indices[0] if len(zero_cross_indices) > 0 else len(ranked_genes) // 2

    if p_adjust >= 0.05:
        return 0
    if direction == 'up' and max_es_idx > zero_cross_idx - len(ranked_genes)/20:
        return 0
    if direction == 'down' and max_es_idx < zero_cross_idx + len(ranked_genes)/20:
        return 0

    # Create GSEA plot
    plt.style.use('ggplot')
    plt.figure(figsize=(8, 8))
    
    x = np.arange(len(ranked_genes))
    plt.plot(x, running_sum, color='blue' if direction == 'up' else 'red', linewidth=2)
    plt.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    plt.axvline(x=max_es_idx, color='blue' if direction == 'up' else 'red', linestyle='--', alpha=0.7)
    plt.axvline(x=zero_cross_idx, color='gray', linestyle='--', alpha=0.7)
    
    # Mark hits
    v_max = np.max(running_sum)
    v_min = np.min(running_sum)
    v_hit = (v_max - v_min) * (1 if np.abs(v_max) > np.abs(v_min) else -1) / 20
    plt.vlines(hit_indices, 0, v_hit, color='black', alpha=0.2)
    
    # Customize
    plt.ylabel('Enrichment Score', fontsize=12)
    plt.xlabel('Rank of Genes', fontsize=12)
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
                        fontsize=12,
                        handlelength=0)

    # Save GSEA plot
    plt.tight_layout()
    gsea_output = os.path.join(plots_dir, f'{direction}_{gene_set_name}_gsea.png')
    plt.savefig(gsea_output, dpi=600, bbox_inches='tight')
    gsea_output = os.path.join(plots_dir, f'{direction}_{gene_set_name}_gsea.pdf')
    plt.savefig(gsea_output, dpi=600, bbox_inches='tight')
    plt.close()
    
    # Create bar plot
    hit_genes = [ranked_genes[i] for i in hit_indices]
    if direction == 'up':
        selected_hits = [gene for i, gene in enumerate(hit_genes) if hit_indices[i] < max_es_idx]
    else:
        selected_hits = [gene for i, gene in enumerate(hit_genes) if hit_indices[i] > max_es_idx]
        selected_hits.reverse()
    
    if selected_hits:
        p_values = [abs(correlations[gene]) for gene in selected_hits]
        plt.figure(figsize=(16, 8))
        bars = plt.bar(range(len(selected_hits)), p_values, 
                      color='blue' if direction == 'up' else 'red', alpha=0.7)
        
        # Customize
        plt.xticks(range(len(selected_hits)), selected_hits, rotation=90)
        plt.ylabel('-log10 p-value', fontsize=12)
        plt.xlabel('Leading Genes', fontsize=12)
        plt.xlim(-1, len(selected_hits))

        # Text box
        gene_set_text = '_'.join(gene_set_name.split('_')[1:])
        gene_cnt_text = f'Gene Count = {len(selected_hits)}'
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, fc='none', ec='none', label=gene_set_text),
            plt.Rectangle((0, 0), 1, 1, fc='none', ec='none', label=gene_cnt_text)
        ]
        legend = plt.legend(handles=legend_elements, 
                            loc='upper right',
                            frameon=False,
                            fontsize=12,
                            handlelength=0)
        
        # Save bar plot
        plt.tight_layout()
        bar_output = os.path.join(plots_dir, f'{direction}_{gene_set_name}_bar.png')
        plt.savefig(bar_output, dpi=600, bbox_inches='tight')
        bar_output = os.path.join(plots_dir, f'{direction}_{gene_set_name}_bar.pdf')
        plt.savefig(bar_output, dpi=600, bbox_inches='tight')
        plt.close()
    
    return 1

def GSEA(summarize_file, group_column, group1, group2, gmt_file):    
    # Read data
    print(f"[Read Data] ...", end='\r')
    df = pd.read_csv(summarize_file)
    start_gene_index = df.columns.get_loc('START_GENE')
    metadata = df.iloc[:, :start_gene_index + 1]
    gene_data = df.iloc[:, start_gene_index + 1:]
    print("[Read Data] Complete    ")

    # Filter and normalize gene
    print(f"[Filter Genes] ...", end='\r')
    expressed_genes = gene_data.columns[gene_data.mean() > 1]
    gene_data = gene_data[expressed_genes]
    high_var_genes = gene_data.columns[gene_data.var() > 0]
    selected_gene_data = gene_data[high_var_genes]
    target_median = selected_gene_data.median(axis=1).median(axis=0)
    scale_factors = target_median / selected_gene_data.median(axis=1)
    expression_data = selected_gene_data.multiply(scale_factors, axis=0)
    print(f"[Filter Genes] {expression_data.shape[1]}    ")

    # Load GeneSets
    print(f"[Load GeneSets] ...", end='\r')
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            gene_set_name = parts[0]
            genes = set(parts[2:])
            gene_sets[gene_set_name] = genes
    print(f"[Load GeneSets] {len(gene_sets)}       ")   

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

    up_results = []
    down_results = []
    total_sets = len(gene_sets)
    for i, (gene_set_name, gene_set) in enumerate(gene_sets.items(), 1):
        es, p_value, direction = calculate_ks_test(ranked_genes, gene_set, 'up')
        up_results.append((gene_set_name, es * direction, p_value))
        es, p_value, direction = calculate_ks_test(ranked_genes, gene_set, 'down')
        down_results.append((gene_set_name, es * direction, p_value))
        print(f"[KS Test] {i}/{total_sets}    ", end='\r')
    
    # Adjust p-values
    _, up_adjusted_pvalues, _, _ = multipletests([r[2] for r in up_results], method='fdr_bh')
    up_final_results = [(r[0], r[1], r[2], adj_p) for r, adj_p in zip(up_results, up_adjusted_pvalues)]
    _, down_adjusted_pvalues, _, _ = multipletests([r[2] for r in down_results], method='fdr_bh')
    down_final_results = [(r[0], r[1], r[2], adj_p) for r, adj_p in zip(down_results, down_adjusted_pvalues)]
    print("[KS Test] Complete                 ")

    # Save results
    print("[Save Results] ...", end='\r')
    geneset_name = os.path.splitext(os.path.basename(gmt_file))[0]
    output_dir = os.path.dirname(summarize_file)

    up_results_df = pd.DataFrame(up_final_results, columns=['gene_set', 'enrichment_score', 'p_value', 'p_adjust'])
    up_results_df = up_results_df.sort_values('p_value')
    up_results_file = os.path.join(output_dir, f'GSEA_{"+".join(group1)}_vs_{"+".join(group2)}_{geneset_name}_up.csv')
    up_results_df.to_csv(up_results_file, index=False)

    down_results_df = pd.DataFrame(down_final_results, columns=['gene_set', 'enrichment_score', 'p_value', 'adjusted_pvalue'])
    down_results_df = down_results_df.sort_values('p_value')
    down_results_file = os.path.join(output_dir, f'GSEA_{"+".join(group1)}_vs_{"+".join(group2)}_{geneset_name}_down.csv')
    down_results_df.to_csv(down_results_file, index=False)
    print("[Save Results] Complete    ")

    print("[Create Plots] ...", end='\r')
    plots_dir = os.path.join(output_dir, f'GSEA_{"+".join(group1)}_vs_{"+".join(group2)}_{geneset_name}_plots')
    os.makedirs(plots_dir, exist_ok=True)

    significant_up = 0
    for _, row in up_results_df.iterrows():
        gene_set_name = row['gene_set']
        p_adjust = row['p_adjust']
        significant_up += plot_gsea_bar(ranked_genes, gene_sets[gene_set_name], gene_set_name, p_adjust, 
                                               correlations, plots_dir, direction='up')
    
    significant_down = 0
    for _, row in down_results_df.iterrows():
        gene_set_name = row['gene_set']
        p_adjust = row['adjusted_pvalue']
        significant_down += plot_gsea_bar(ranked_genes, gene_sets[gene_set_name], gene_set_name, p_adjust, 
                                                 correlations, plots_dir, direction='down')
    
    print(f"[Create Plots] Up: {significant_up} / Down: {significant_down}    ")


if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("Usage: python3 GSEA.py <cohort/summarize.csv> <group_column> <group_a1> <group_a2> ... -- <group_b1> <group_b2> ... <geneset.gmt>")
        sys.exit(1)
    
    summarize_file = sys.argv[1]
    group_column = sys.argv[2]
    separator_index = sys.argv.index('--')
    group1 = sys.argv[3:separator_index]
    group2 = sys.argv[separator_index+1:-1]
    gmt_file = sys.argv[-1]
    
    GSEA(summarize_file, group_column, group1, group2, gmt_file)