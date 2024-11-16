import pandas as pd

def process_csv(input_file, output_file):
    # Read the TSV file, skipping the first 2 metadata rows
    df = pd.read_csv(input_file, sep='\t', skiprows=[1, 2])
    
    # Get column names
    columns = df.columns.tolist()
    
    # Swap first two columns
    columns[0], columns[1] = columns[1], columns[0]
    df = df[columns]
    
    # Replace hyphens with underscores in the TCGA identifiers
    df.iloc[:, 0] = df.iloc[:, 0].apply(lambda x: x.replace('-', '_') if isinstance(x, str) and x.startswith('TCGA-') else x)
    
    # Save to new file
    df.to_csv(output_file, sep='\t', index=False)

# Example usage
input_file = 'LIHC.tsv'  # Replace with your input file path
output_file = 'LIHC_fixed.tsv'  # Replace with desired output file path

process_csv(input_file, output_file)