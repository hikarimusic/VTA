import os
import sys
import pandas as pd
import glob
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

def process_profiles(directory, cohort_file):
    output_file = os.path.join(directory, 'profile.csv')
   
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

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 profile.py <directory> <cohort>")
        sys.exit(1)
   
    directory = sys.argv[1]
    cohort_file = sys.argv[2]
    
    if not os.path.exists(cohort_file):
        print(f"Error: {cohort_file} not found")
        sys.exit(1)
    
    profile_df = process_profiles(directory, cohort_file)
    ddx11l2_index = profile_df.columns.get_loc('DDX11L2')
    for column in profile_df.columns[1:ddx11l2_index]:
        if profile_df[column].dtype == 'object' or pd.api.types.is_categorical_dtype(profile_df[column]):
            pca_output_file = os.path.join(directory, f'PCA_{column}.png')
            generate_pca_plot(profile_df, pca_output_file, column)