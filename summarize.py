import os
import sys
import pandas as pd
import numpy as np

def summarize(cohort_file, profile_dir, value_type):
    # Process TSV
    print(f"[Process TSV] ...    ", end='\r')
    file_extension = os.path.splitext(cohort_file)[1].lower()
    separator = '\t' if file_extension == '.tsv' else ','
    cohort_df = pd.read_csv(cohort_file,sep=separator)
    sample_ids = cohort_df.iloc[:, 0].tolist()
   
    gene_value_dict = {}
    gene_order = []
    for i, sample_id in enumerate(sample_ids, 1):
        tsv_file = os.path.join(profile_dir, f"{sample_id}.tsv")
        if not os.path.exists(tsv_file):
            print(f"[Warning] {sample_id}.tsv not found                 ")
            continue
       
        print(f"[Process TSV] {i}/{len(sample_ids)}                 ", end='\r')
        df = pd.read_csv(tsv_file, sep='\t', comment='#')
        df = df[~df['gene_id'].str.startswith('N_')]
        gene_value_dict[sample_id] = dict(zip(df['gene_name'], df[value_type]))
        if not gene_order:
            gene_order = df['gene_name'].tolist()
    gene_value_df = pd.DataFrame({gene: [gene_value_dict[sample].get(gene, 0) for sample in sample_ids if sample in gene_value_dict]
                                  for gene in gene_order})
    print("[Process TSV] Complete")

    # Save result
    print("[Save Result] ...    ", end='\r')
    output_dir = os.path.splitext(cohort_file)[0]
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'summary.csv')
    output_file2 = os.path.join(output_dir, 'cohort.csv')
   
    # Add START_GENE column
    cohort_df = cohort_df[cohort_df.iloc[:, 0].isin(gene_value_dict.keys())]
    cohort_df = cohort_df.reset_index(drop=True)
    cohort_df['START_GENE'] = ""
    result_df = pd.concat([cohort_df, gene_value_df], axis=1)
   
    result_df.to_csv(output_file, index=False)
    cohort_df.to_csv(output_file2, index=False)
    print("[Save Result] Complete                 ")
    

if __name__ == "__main__":
    '''Command: python3 summarize.py <cohort.csv> <profile_dir/> <value_type>'''
    cohort_file = sys.argv[1]
    profile_dir = sys.argv[2]
    value_type = sys.argv[3]
   
    summarize(cohort_file, profile_dir, value_type)