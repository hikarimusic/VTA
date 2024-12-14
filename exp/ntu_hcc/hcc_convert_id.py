import pandas as pd

def load_gsm_srr_mapping(mapping_file):
    """Load the GSM to SRR mapping from the list file."""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            gsm, srr = line.strip().split()
            mapping[gsm] = srr
    return mapping

def process_hcc_data(family_file, mapping_file):
    """Process HCC family data with GSM to SRR conversion and column renaming."""
    # Load the GSM to SRR mapping
    gsm_to_srr = load_gsm_srr_mapping(mapping_file)
    
    # Read the family CSV file
    df = pd.read_csv(family_file)
    
    # Convert GSM IDs to SRR IDs
    df['sample_id'] = df['sample_id'].map(gsm_to_srr)
    
    # Define the new column names
    column_mapping = {
        'sample_id': 'Sample',
        'etiology': 'Etiology',
        'ctnnb1 status': 'CTNNB1',
        'efs days': 'EFS_Days',
        'efs event': 'EFS_event',
        'os days': 'OS_days',
        'os event': 'OS_event'
    }
    
    # Rename the columns
    df = df.rename(columns=column_mapping)
    
    return df

def main():
    # Define input files
    family_file = 'hcc_family.csv'
    mapping_file = 'hcc_list.txt'
    output_file = 'hcc_cohort.csv'
    
    # Process the data
    result_df = process_hcc_data(family_file, mapping_file)
    
    # Save to CSV
    result_df.to_csv(output_file, index=False)
    print(f"Processed data saved to {output_file}")

if __name__ == "__main__":
    main()