import pandas as pd
import numpy as np

# Read the input TSV file
df = pd.read_csv('clinical.tsv', sep='\t')

# Drop duplicates based on case_id, keeping the first occurrence
df = df.drop_duplicates(subset=['case_id'], keep='first')

# Create a copy of the filtered dataframe to avoid warnings
df = df[df['vital_status'].isin(['Dead', 'Alive'])].copy()

# Create new 'event' column (0 for Alive, 1 for Dead)
df.loc[:, 'event'] = (df['vital_status'] == 'Dead').astype(int)

# Create new 'time' column based on conditions
df.loc[:, 'time'] = np.where(
    df['event'] == 1,
    df['days_to_death'],
    df['days_to_last_follow_up']
)

# Convert '--' to NaN in the time column if present
df.loc[:, 'time'] = df['time'].replace("'--", np.nan)

# Convert time to numeric, coercing errors to NaN
df.loc[:, 'time'] = pd.to_numeric(df['time'], errors='coerce')

# Get the first two column names
cols = df.columns.tolist()
first_col = cols[0]
second_col = cols[1]

# Swap the first two columns
cols[0], cols[1] = cols[1], cols[0]

# Reorder the DataFrame columns
df = df[cols]

# Convert hyphens to underscores in the first column (case_submitter_id)
df.iloc[:, 0] = df.iloc[:, 0].str.replace('-', '_')

# Save the processed data to the output TSV file
df.to_csv('clinical_fine.tsv', sep='\t', index=False)

# Print statistics
original_rows = len(pd.read_csv('clinical.tsv', sep='\t'))
print(f"Original number of rows: {original_rows}")
print(f"Number of rows after deduplication: {len(df[df['vital_status'].isin(['Dead', 'Alive'])])}")
print(f"Number of rows in final output: {len(df)}")
print("\nSample of event and time columns:")
print(df[['vital_status', 'event', 'time']].head())
print("\nSample of converted case submitter IDs:")
print(df.iloc[:5, 0])  # Show first 5 entries of first column
print("\nColumns in output file:")
print(df.columns.tolist()[:5], "...")  # Show first 5 columns