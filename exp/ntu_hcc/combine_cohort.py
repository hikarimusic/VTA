import pandas as pd
import numpy as np

# Read LIVER cohort
liver_df = pd.read_csv('LIVER_cohort.csv')
liver_df['Cohort'] = 'LIVER'

# Read NTUHCC cohort
ntuhcc_df = pd.read_csv('NTUHCC_cohort.csv')
ntuhcc_df['Cohort'] = 'NTUHCC'

# Read LIHC cohort
lihc_df = pd.read_csv('LIHC_cohort.tsv', sep='\t')
lihc_df['Cohort'] = 'LIHC'
lihc_df['Sample'] = lihc_df['bcr_patient_barcode']

# Fill missing values with 'none'
liver_df = liver_df.fillna('none')
ntuhcc_df = ntuhcc_df.fillna('none')
lihc_df = lihc_df.fillna('none')

# Combine all dataframes
# First, let's identify common columns if any
all_dfs = [liver_df, ntuhcc_df, lihc_df]
combined_df = pd.concat(all_dfs, ignore_index=True, sort=False)

# Ensure all columns have 'none' for missing values
combined_df = combined_df.fillna('none')

# Save the combined dataset
combined_df.to_csv('combined_cohorts.csv', index=False)

# Print some information about the combined dataset
print(f"Total number of samples: {len(combined_df)}")
print("\nNumber of samples per cohort:")
print(combined_df['Cohort'].value_counts())
print("\nColumns in combined dataset:")
print(combined_df.columns.tolist())