import pandas as pd
import sys

def csv_to_gct(input_csv, output_gct):
    """
    Convert a CSV file containing gene expression data to GCT format.
    Gene columns start after the 'START_GENE' column.
    Sample names are in the first column.
    
    Parameters:
    input_csv (str): Path to input CSV file
    output_gct (str): Path to output GCT file
    """
    try:
        # Read the CSV file
        df = pd.read_csv(input_csv, low_memory=False)
        
        # Get sample names from the first column
        sample_names = df.iloc[:, 0].values
        
        # Find the index where gene columns start
        start_gene_index = df.columns.get_loc('START_GENE')
        
        # Get gene columns (all columns after START_GENE)
        gene_cols = df.columns[start_gene_index + 1:]
        
        # Prepare data for GCT format
        # Get gene expression data
        gene_data = df[gene_cols]
        
        # Transpose the data and set sample names as columns
        gct_data = gene_data.T
        gct_data.columns = sample_names
        
        # Add Description column (using 'NA' as placeholder)
        gct_data.insert(0, 'Description', 'NA')
        
        # Use gene names as index/Name column
        gct_data.index.name = 'Name'
        
        # Get dimensions for the second line of GCT file
        num_rows, num_cols = gct_data.shape
        num_cols = num_cols - 1  # Subtract 1 to exclude Description column
        
        # Write to GCT file
        with open(output_gct, 'w') as f:
            # Write version
            f.write("#1.2\n")
            
            # Write dimensions
            f.write(f"{num_rows}\t{num_cols}\n")
            
            # Convert DataFrame to GCT format and write
            gct_data.to_csv(f, sep='\t', mode='a')
            
        print(f"Successfully converted {input_csv} to GCT format. Output saved as {output_gct}")
        print(f"Number of genes: {num_rows}")
        print(f"Number of samples: {num_cols}")
        
    except KeyError:
        print("Error: 'START_GENE' column not found in the CSV file")
        sys.exit(1)
    except Exception as e:
        print(f"Error occurred during conversion: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    # Example usage
    input_file = "summary_cluster.csv"
    output_file = "summary_cluster.gct"
    csv_to_gct(input_file, output_file)