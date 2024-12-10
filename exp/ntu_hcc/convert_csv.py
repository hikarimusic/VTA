import pandas as pd
import sys

def gmt_to_csv(input_gmt, output_csv, original_csv):
    """
    Convert a GMT file to CSV format, replacing gene columns with pathway enrichment scores.
    
    Parameters:
    input_gmt (str): Path to input GMT file
    output_csv (str): Path to output CSV file
    original_csv (str): Path to original CSV file to get the structure
    """
    try:
        # Read the original CSV to get the structure
        original_df = pd.read_csv(original_csv, low_memory=False)
        start_gene_index = original_df.columns.get_loc('START_GENE')
        metadata = original_df.iloc[:, :start_gene_index + 1]
        
        # Read the GMT file
        with open(input_gmt, 'r') as f:
            # Skip the first two lines (version and dimensions)
            next(f)
            next(f)
            # Read the rest as a CSV
            gmt_df = pd.read_csv(f, sep='\t')
        
        # Process pathway names
        gmt_df['Pathway'] = gmt_df['Name'].apply(lambda x: '_'.join(x.split('_')[1:]))
        
        # Remove 'Name' and 'Description' columns and set Pathway as index
        score_df = gmt_df.drop(['Name', 'Description'], axis=1)
        score_df.index = gmt_df['Pathway']
        
        # Create the new DataFrame with metadata
        final_df = metadata.copy()
        
        # For each sample in metadata
        for idx, row in final_df.iterrows():
            sample_id = row.iloc[0]  # Get sample ID from first column
            if sample_id in score_df.columns:
                # Add pathway scores for this sample
                for pathway in score_df.index:
                    final_df.loc[idx, pathway] = score_df.loc[pathway, sample_id]
        
        # Save to CSV
        final_df.to_csv(output_csv, index=False)
        
        print(f"Successfully converted {input_gmt} to CSV format.")
        print(f"Output saved as {output_csv}")
        print(f"Number of pathways added: {len(score_df.index)}")
        print(f"Number of samples processed: {len(metadata)}")
        
    except Exception as e:
        print(f"Error occurred during conversion: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    # Example usage
    input_file = "summary_ssgsea.gmt"
    output_file = "summary_ssgsea.csv"
    original_file = "summary_cluster.csv"
    gmt_to_csv(input_file, output_file, original_file)