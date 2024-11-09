import sys
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

def generate_pca_plot(summarize_file, group_columns):
    print(f"[Reading Data] ...", end='\r')
    if not os.path.exists(summarize_file):
        print(f"\n[PCA plot] Error: {summarize_file} not found")
        sys.exit(1)
   
    df = pd.read_csv(summarize_file)
    start_gene_index = df.columns.get_loc('START_GENE')
    features = df.iloc[:, start_gene_index + 1:]
    labels = df[group_columns]
    print("[Reading Data] Complete")
   
    # Normalize the features
    print(f"[PCA plot] ...", end='', flush=True)
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)
   
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_features)
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df = pd.concat([pca_df, labels], axis=1)
   
    # Create plots
    mpl.style.use('ggplot')
    for group_column in group_columns:
        plt.figure(figsize=(10, 10))
        sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=group_column, palette='deep', s=100)
       
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance explained)', fontsize=14)
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance explained)', fontsize=14)
        plt.legend(title=group_column, title_fontsize='13', fontsize='11')
       
        output_dir = os.path.dirname(summarize_file)
        output_file = os.path.join(output_dir, f'PCA_{group_column.replace(" ", "_")}.pdf')
        plt.savefig(output_file, format='pdf', dpi=600, bbox_inches='tight')
        output_file = os.path.join(output_dir, f'PCA_{group_column.replace(" ", "_")}.png')
        plt.savefig(output_file, format='png', dpi=600, bbox_inches='tight')
        plt.close()
   
    print(f"\r[PCA plot] Complete                   ", flush=True)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 PCA.py <cohort/summarize.csv> <group_column1> <group_column2> ...")
        sys.exit(1)
   
    summarize_file = sys.argv[1]
    group_columns = sys.argv[2:]
   
    generate_pca_plot(summarize_file, group_columns)