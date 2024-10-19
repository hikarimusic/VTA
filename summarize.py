import os
import sys
import pandas as pd
import numpy as np

def process_profiles(cohort_file, profile_dir, value_type, start_gene):
    output_dir = os.path.splitext(cohort_file)[0]
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'summary.csv')
    output_file2 = os.path.join(output_dir, 'cohort.csv')
   
    cohort_df = pd.read_csv(cohort_file)
    sample_ids = cohort_df.iloc[:, 0].tolist()
   
    gene_value_dict = {}
    gene_order = []
    for i, sample_id in enumerate(sample_ids, 1):
        tsv_file = os.path.join(profile_dir, f"{sample_id}.tsv")
        if not os.path.exists(tsv_file):
            print(f"[Warning] {sample_id}.tsv not found")
            continue
       
        print(f"[Process TSV] {i}/{len(sample_ids)}", end='\r')
        df = pd.read_csv(tsv_file, sep='\t')
        gene_value_dict[sample_id] = dict(zip(df['gene_name'], df[value_type]))
        if not gene_order:
            gene_order = df['gene_name'].tolist()
    print("[Process TSV] Complete")
   
    print("[Create Summary] ...", end='\r')
    gene_value_df = pd.DataFrame({gene: [gene_value_dict[sample].get(gene, 0) for sample in sample_ids if sample in gene_value_dict]
                                  for gene in gene_order})
    print("[Create Summary] Complete")
   
    result_df = pd.concat([cohort_df[cohort_df.iloc[:, 0].isin(gene_value_dict.keys())], gene_value_df], axis=1)
   
    # Add START_GENE column
    start_gene_index = result_df.columns.get_loc(start_gene)
    result_df.insert(start_gene_index, 'START_GENE', "")
   
    print("[Save Result] ...", end='\r')
    result_df.to_csv(output_file, index=False)
    cohort_df.to_csv(output_file2, index=False)
    print("[Save Result] Complete")
   
    return result_df

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 summarize.py <cohort.csv> <profile_dir/> <start_gene> <value_type>")
        sys.exit(1)
   
    cohort_file = sys.argv[1]
    profile_dir = sys.argv[2]
    start_gene = sys.argv[3]
    value_type = sys.argv[4]
   
    if not os.path.exists(cohort_file):
        print(f"Error: {cohort_file} not found")
        sys.exit(1)
   
    if not os.path.exists(profile_dir):
        print(f"Error: {profile_dir} not found")
        sys.exit(1)
   
    process_profiles(cohort_file, profile_dir, value_type, start_gene)