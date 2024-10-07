import pandas as pd

# Read the file
df = pd.read_csv('GSE141198_TLCN.subset1.readcount.txt', sep='\t')

# Transpose the dataframe
df_transposed = df.set_index('GeneID').T

# Reset the index to turn sample names into a column
df_transposed = df_transposed.reset_index()

# Rename the 'index' column to 'sample'
df_transposed = df_transposed.rename(columns={'index': 'sample'})

# Save the result to a new CSV file
df_transposed.to_csv('hcc_converted.csv', index=False)

print("Conversion complete. Output saved to 'hcc_converted.csv'")