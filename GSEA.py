import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def calculate_correlation(group1_data, group2_data):
    if group1_data.empty or group2_data.empty:
        return 0
    

    combined_data = np.concatenate([group1_data, group2_data])
    labels = np.concatenate([np.zeros(len(group1_data)), np.ones(len(group2_data))])
    correlation, p_value = stats.pearsonr(combined_data, labels)
    
    return correlation
    
def calculate_correlation(group1_data, group2_data):
    if group1_data.empty or group2_data.empty:
        return 0

    # mean1 = group1_data.mean()
    # mean2 = group2_data.mean()
    # log2_fold_change = np.log2((abs(mean2) + 1e-8) / (abs(mean1) + 1e-8))
    # _, p_value = stats.ttest_ind(group1_data, group2_data)

    # return log2_fold_change * -np.log10(p_value)

    combined_data = np.concatenate([group1_data, group2_data])
    labels = np.concatenate([np.zeros(len(group1_data)), np.ones(len(group2_data))])
    correlation, p_value = stats.pearsonr(combined_data, labels)

    return correlation

def calculate_es(ranked_genes, gene_set, weights):
    gene_set = [gene for gene in gene_set if gene in weights]
    N = len(ranked_genes)
    Nh = len(gene_set)
    
    sum_weights = sum(abs(weights[gene]) for gene in gene_set)
    running_sum = 0
    max_es = 0
    for i in range(N):
        gene = ranked_genes[i]
        if gene in gene_set:
            running_sum += (abs(weights[gene]) / sum_weights)
            # running_sum += 1 / Nh
        else:
            running_sum -= 1 / (N - Nh)
        
        if abs(running_sum) > abs(max_es):
            max_es = running_sum

    return max_es

def main(summarize_file, group_column, group1, group2, gmt_file): 
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

    # Load GeneSets
    print(f"[Loading GeneSets] ...", end='\r')
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            gene_set_name = parts[0]
            genes = set(parts[2:])
            gene_sets[gene_set_name] = genes
    print(f"[Loading GeneSets] {len(gene_sets)}       ")   

    # Perform GSEA
    print(f"[GSEA Start] ...", end='\r')
    phenotypes = metadata[group_column].map(lambda x: 1 if x in group1 else (-1 if x in group2 else 0))
    expression_data = expression_data.loc[phenotypes != 0]
    phenotypes = phenotypes[phenotypes != 0]
    genes = expression_data.columns.tolist()

    group1_data = expression_data[phenotypes == 1]
    group2_data = expression_data[phenotypes == -1]
    correlations = {gene: calculate_correlation(group1_data[gene], group2_data[gene]) for gene in genes}
    ranked_genes = sorted(correlations, key=correlations.get, reverse=True)
    weights = {gene: correlations[gene] for gene in ranked_genes}
    
    real_es = {gene_set_name: calculate_es(ranked_genes, gene_set, weights) 
               for gene_set_name, gene_set in gene_sets.items()}
    
    # Permutation test
    num_permutations = 1000
    permutations = [np.random.permutation(phenotypes) for _ in range(num_permutations)]
    null_es = {gene_set_name: [] for gene_set_name in gene_sets}
    
    for i, permuted_phenotypes in enumerate(permutations, 1):
        print(f"[Permutation] {i}/{num_permutations}  ", end='\r')
        permuted_group1_data = expression_data[permuted_phenotypes == 1]
        permuted_group2_data = expression_data[permuted_phenotypes == -1]
        permuted_correlations = {gene: calculate_correlation(permuted_group1_data[gene], permuted_group2_data[gene]) 
                                 for gene in genes}
        permuted_ranked_genes = sorted(permuted_correlations, key=permuted_correlations.get, reverse=True)
        permuted_weights = {gene: permuted_correlations[gene] for gene in permuted_ranked_genes}
        
        for gene_set_name, gene_set in gene_sets.items():
            es = calculate_es(permuted_ranked_genes, gene_set, permuted_weights)
            null_es[gene_set_name].append(es)
    
    print("[Permutation] Complete      ")

    # Calculate p-values
    print("[p-values] ...", end='\r')
    results = []
    for gene_set_name, es in real_es.items():
        null_distribution = null_es[gene_set_name]
        if es >= 0:
            p_value = sum(null_es >= es for null_es in null_distribution) / num_permutations
        else:
            p_value = sum(null_es <= es for null_es in null_distribution) / num_permutations
        results.append((gene_set_name, es, p_value))

    _, adjusted_pvalues, _, _ = multipletests([r[2] for r in results], method='fdr_bh')
    final_results = [(r[0], r[1], r[2], adj_p) for r, adj_p in zip(results, adjusted_pvalues)]
    print("[p-values] Complete")

    # Save results
    print("[Saving Results] ...", end='\r')
    results_df = pd.DataFrame(final_results, columns=['gene_set', 'ES', 'p_value', 'adjusted_pvalue'])
    results_df = results_df.sort_values('p_value')
    
    output_dir = os.path.dirname(summarize_file)
    geneset_name = os.path.splitext(os.path.basename(gmt_file))[0]
    results_file = os.path.join(output_dir, f'GSEA_{"+".join(group1)}_vs_{"+".join(group2)}_{geneset_name}.csv')
    results_df.to_csv(results_file, index=False)
    print("[Saving Results] Complete")

if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("Usage: python3 GSEA.py <cohort/summarize.csv> <group_column> <group1a> <group1b> ... -- <group2a> <group2b> ... <geneset.gmt>")
        sys.exit(1)
    
    summarize_file = sys.argv[1]
    group_column = sys.argv[2]
    
    try:
        separator_index = sys.argv.index('--')
        group1 = sys.argv[3:separator_index]
        group2 = sys.argv[separator_index+1:-1]
        gmt_file = sys.argv[-1]
    except ValueError:
        print("Error: Missing '--' separator between groups")
        sys.exit(1)
    
    if not group1 or not group2:
        print("Error: Both groups must contain at least one class")
        sys.exit(1)
    
    if not os.path.exists(gmt_file):
        print(f"Error: {gmt_file} not found")
        sys.exit(1)
    
    main(summarize_file, group_column, group1, group2, gmt_file)