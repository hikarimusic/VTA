import os
import sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster import hierarchy
from statsmodels.stats.multitest import multipletests

def process_profiles(directory, cohort_file, output_dir):
    output_file = os.path.join(output_dir, 'profile.csv')
   
    if os.path.exists(output_file):
        print(f"[Skip] Profile already exists")
        return pd.read_csv(output_file)
   
    cohort_df = pd.read_csv(cohort_file)
    sample_ids = cohort_df.iloc[:, 0].tolist()
   
    gene_tpm_dict = {}
    gene_order = []
    for i, sample_id in enumerate(sample_ids, 1):
        tsv_file = os.path.join(directory, f"{sample_id}.tsv")
        if not os.path.exists(tsv_file):
            print(f"\r[Warning] {sample_id}.tsv not found")
            continue
       
        print(f"\r[Process TSV] {i}/{len(sample_ids)}                    ", end='', flush=True)
        df = pd.read_csv(tsv_file, sep='\t')
        gene_tpm_dict[sample_id] = dict(zip(df['gene_name'], df['tpm']))
        if not gene_order:
            gene_order = df['gene_name'].tolist()
    print("\r[Process TSV] Complete                    ", flush=True)
   
    print("\r[Create TPM] ...                    ", end='', flush=True)
    gene_tpm_df = pd.DataFrame({gene: [gene_tpm_dict[sample].get(gene, 0) for sample in sample_ids if sample in gene_tpm_dict]
                                for gene in gene_order})
    print("\r[Create TPM] Complete                    ", flush=True)
   
    result_df = pd.concat([cohort_df[cohort_df.iloc[:, 0].isin(gene_tpm_dict.keys())], gene_tpm_df], axis=1)
   
    print("\r[Save Result] ...                    ", end='', flush=True)
    result_df.to_csv(output_file, index=False)
    print("\r[Save Result] Complete                    ", flush=True)
   
    return result_df

def generate_pca_plot(df, output_file, hue_column):
    if os.path.exists(output_file):
        print(f"[Skip] PCA already exists")
        return
    
    print(f"\r[Generate PCA] {hue_column} ...                    ", end='', flush=True)
    ddx11l2_index = df.columns.get_loc('DDX11L2')
    features = df.iloc[:, ddx11l2_index:]
    labels = df[hue_column]
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_features)
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df[hue_column] = labels
    mpl.style.use('ggplot')
    plt.figure(figsize=(12, 8))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=hue_column, palette='deep', s=100)
    plt.title(f'PCA of RNA Expression Profiles (Colored by {hue_column})', fontsize=20)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance explained)', fontsize=14)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance explained)', fontsize=14)
    plt.legend(title=hue_column, title_fontsize='13', fontsize='11')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
   
    print(f"\r[Generate PCA] {hue_column} Complete                    ", flush=True)

def read_comparison_config(config_file):
    comparisons = []
    current_comparison = {}
    
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Column:'):
                if current_comparison:
                    comparisons.append(current_comparison)
                    current_comparison = {}
                current_comparison['column'] = line.split(':', 1)[1].strip()
            elif line.startswith('Group1:'):
                current_comparison['group1'] = [g.strip() for g in line.split(':', 1)[1].split(',')]
            elif line.startswith('Group2:'):
                current_comparison['group2'] = [g.strip() for g in line.split(':', 1)[1].split(',')]
    
    if current_comparison:
        comparisons.append(current_comparison)
    
    return comparisons

def perform_differential_expression(df, group_column, group1, group2, output_dir):
    print(f"\r[Differential Expression] Analyzing {group_column}: {group1} vs {group2}...                    ", end='', flush=True)
    
    # Separate metadata and expression data
    metadata = df.iloc[:, :df.columns.get_loc('DDX11L2')]
    expression_data = df.iloc[:, df.columns.get_loc('DDX11L2'):]
    
    # Filter data for the specified groups
    group1_data = expression_data[metadata[group_column].isin(group1)]
    group2_data = expression_data[metadata[group_column].isin(group2)]
    
    # print(f"\nGroup 1 ({'+'.join(group1)}) sample count: {len(group1_data)}")
    # print(f"Group 2 ({'+'.join(group2)}) sample count: {len(group2_data)}")
    
    if len(group1_data) == 0 or len(group2_data) == 0:
        print(f"Error: One or both groups have no samples. Skipping differential expression analysis.")
        return None
    
    # Perform t-test for each gene
    pvalues = []
    log2_fold_changes = []
    genes = []
    
    for gene in expression_data.columns:
        # Check if all values are the same in both groups
        if group1_data[gene].nunique() == 1 and group2_data[gene].nunique() == 1 and group1_data[gene].iloc[0] == group2_data[gene].iloc[0]:
            continue
        
        t_stat, p_value = stats.ttest_ind(group1_data[gene], group2_data[gene])
        
        # Handle NaN p-values
        if np.isnan(p_value):
            print(f"Skipping gene {gene}: NaN p-value")
            continue
        
        pvalues.append(p_value)
        genes.append(gene)
        
        # Handle potential division by zero or log of zero
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
    
    total_genes = len(genes)
    # print(f"Total genes analyzed: {total_genes}")
    # print(f"Number of genes skipped: {len(expression_data.columns) - total_genes}")
    
    # Bonferroni correction
    bonferroni_adjusted = np.minimum(np.array(pvalues) * total_genes, 1.0)
    
    
    # Create DataFrame with results
    results = pd.DataFrame({
        'gene': genes,
        'log2_fold_change': log2_fold_changes,
        'p_value': pvalues,
        'bonferroni_adjusted': bonferroni_adjusted,
    })
    
    # Sort by Bonferroni adjusted p-value
    results = results.sort_values('bonferroni_adjusted')
    
    # Save results
    output_file = os.path.join(output_dir, f'differential_expression_{group_column}_{"+".join(group1)}_vs_{"+".join(group2)}.csv')
    results.to_csv(output_file, index=False)
    
    # Generate volcano plot
    plt.style.use('ggplot')
    plt.figure(figsize=(10, 8))
    
    # Set thresholds
    log2_fc_threshold = 1  # log2 fold change threshold
    p_value_threshold = 0.05  # p-value threshold
    
    # Create a color column
    results['color'] = 'grey'
    results.loc[(results['log2_fold_change'] > log2_fc_threshold) & (results['bonferroni_adjusted'] < p_value_threshold), 'color'] = 'red'
    results.loc[(results['log2_fold_change'] < -log2_fc_threshold) & (results['bonferroni_adjusted'] < p_value_threshold), 'color'] = 'blue'
    
    # Plot points
    plt.scatter(results['log2_fold_change'], -np.log10(results['bonferroni_adjusted']), c=results['color'], alpha=0.5)
    
    # Add threshold lines
    plt.axvline(x=log2_fc_threshold, color='gray', linestyle='--')
    plt.axvline(x=-log2_fc_threshold, color='gray', linestyle='--')
    plt.axhline(y=-np.log10(p_value_threshold), color='gray', linestyle='--')
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Bonferroni Adjusted P-value')
    plt.title(f'Volcano Plot: {group_column} - {"+".join(group2)} vs {"+".join(group1)}')
    
    plot_file = os.path.join(output_dir, f'volcano_plot_{group_column}_{"+".join(group1)}_vs_{"+".join(group2)}.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\r[Differential Expression] Complete: {group_column}: {group1} vs {group2}                    ")
    return results

# def generate_heatmap(df, group_column, groups, n_top_genes=50, output_dir='.'):
#     print(f"\r[Generate Heatmap] Creating for {group_column}: {groups}...                    ", end='', flush=True)
    
#     # Separate metadata and expression data
#     metadata = df.iloc[:, :df.columns.get_loc('DDX11L2')]
#     expression_data = df.iloc[:, df.columns.get_loc('DDX11L2'):]
    
#     # Filter data for the specified groups
#     filtered_metadata = metadata[metadata[group_column].isin(groups)]
#     filtered_expression_data = expression_data.loc[filtered_metadata.index]
    
#     print(f"\nSample count: {len(filtered_expression_data)}")
#     print(f"Group counts: {filtered_metadata[group_column].value_counts().to_dict()}")
    
#     if len(filtered_expression_data) == 0:
#         print(f"Error: No samples found for the specified groups. Skipping heatmap generation.")
#         return
    
#     # Select top N genes based on variance
#     gene_variance = filtered_expression_data.var().sort_values(ascending=False)
#     top_genes = gene_variance.head(n_top_genes).index
    
#     # Prepare data for clustering
#     data_for_clustering = filtered_expression_data[top_genes].T
    
#     # Perform hierarchical clustering
#     linkage = hierarchy.linkage(data_for_clustering, method='average')
    
#     # Create heatmap
#     plt.figure(figsize=(12, 10))
    
#     # Plot dendrogram
#     try:
#         dendrogram = hierarchy.dendrogram(linkage, labels=filtered_metadata[group_column], leaf_rotation=90, leaf_font_size=8)
#     except ValueError as e:
#         print(f"Error: {str(e)}. Skipping dendrogram.")
#         dendrogram = {'leaves': range(len(filtered_metadata))}
    
#     # Reorder data based on dendrogram
#     reordered_data = data_for_clustering.iloc[:, dendrogram['leaves']]
    
#     # Plot heatmap
#     sns.heatmap(reordered_data, cmap='RdBu_r', center=0, xticklabels=False)
    
#     plt.title(f'Gene Expression Heatmap: {group_column} - {", ".join(groups)}')
#     plt.ylabel('Genes')
#     plt.tight_layout()
#     plot_file = os.path.join(output_dir, f'expression_heatmap_{group_column}_{"+".join(groups)}.png')
#     plt.savefig(plot_file, dpi=300, bbox_inches='tight')
#     plt.close()
    
#     print(f"\r[Generate Heatmap] Complete for {group_column}: {groups}                    ")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 profile.py <directory> <cohort> <comparison_config>")
        sys.exit(1)
    
    directory = sys.argv[1]
    cohort_file = sys.argv[2]
    comparison_config = sys.argv[3]
    
    if not os.path.exists(cohort_file):
        print(f"Error: {cohort_file} not found")
        sys.exit(1)
    
    if not os.path.exists(comparison_config):
        print(f"Error: {comparison_config} not found")
        sys.exit(1)
    
    # Create output directory
    cohort_name = os.path.splitext(os.path.basename(cohort_file))[0]
    output_dir = os.path.join(directory, cohort_name)
    os.makedirs(output_dir, exist_ok=True)
    
    profile_df = process_profiles(directory, cohort_file, output_dir)
    
    # Read comparison configuration
    comparisons = read_comparison_config(comparison_config)
    
    for comparison in comparisons:
        column = comparison['column']
        group1 = comparison['group1']
        group2 = comparison['group2']
        
        # Perform differential expression analysis
        diff_expr_results = perform_differential_expression(profile_df, column, group1, group2, output_dir)
        
        # # Generate heatmap for each group
        # generate_heatmap(profile_df, column, group1, n_top_genes=50, output_dir=output_dir)
        # generate_heatmap(profile_df, column, group2, n_top_genes=50, output_dir=output_dir)
        
        # # Generate combined heatmap
        # generate_heatmap(profile_df, column, group1 + group2, n_top_genes=50, output_dir=output_dir)