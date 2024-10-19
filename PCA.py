import sys
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

def generate_pca_plot(cohort_file, group_columns):
    print(f"[Reading Data] ...", end='\r')
    # Construct the path to the summarize.csv file
    cohort_dir = os.path.splitext(cohort_file)[0]
    summarize_file = os.path.join(cohort_dir, 'summarize.csv')
   
    # Check if the summarize.csv file exists
    if not os.path.exists(summarize_file):
        print(f"\n[PCA plot] Error: {summarize_file} not found")
        sys.exit(1)
   
    # Read the summarize.csv file
    df = pd.read_csv(summarize_file)
   
    # Find the index of the START_GENE column
    start_gene_index = df.columns.get_loc('START_GENE')
   
    # Extract features (gene expression data) and labels
    features = df.iloc[:, start_gene_index + 1:]  # +1 to start from the column after START_GENE
    labels = df[group_columns]
    print("[Reading Data] Complete")
   
    # Normalize the features
    print(f"[PCA plot] ...", end='', flush=True)
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)
   
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_features)
   
    # Create a DataFrame with PCA results
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df = pd.concat([pca_df, labels], axis=1)
   
    # Set the style to ggplot
    mpl.style.use('ggplot')
   
    # Create a plot for each group column
    for group_column in group_columns:
        plt.figure(figsize=(10, 10))
        sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=group_column, palette='deep', s=100)
       
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance explained)', fontsize=14)
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance explained)', fontsize=14)
        plt.legend(title=group_column, title_fontsize='13', fontsize='11')
       
        # Save the plot in the same directory as summarize.csv
        output_file = os.path.join(cohort_dir, f'PCA_{group_column.replace(" ", "_")}.pdf')
        plt.savefig(output_file, format='pdf', dpi=600, bbox_inches='tight')
        output_file = os.path.join(cohort_dir, f'PCA_{group_column.replace(" ", "_")}.png')
        plt.savefig(output_file, format='png', dpi=600, bbox_inches='tight')
        plt.close()
   
    print(f"\r[PCA plot] Complete                   ", flush=True)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 PCA.py <cohort_file.csv> <group_column1> <group_column2> ...")
        sys.exit(1)
   
    cohort_file = sys.argv[1]
    group_columns = sys.argv[2:]
   
    generate_pca_plot(cohort_file, group_columns)