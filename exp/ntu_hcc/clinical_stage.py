import pandas as pd

def standardize_pathologic_stages(df):
    """
    Standardizes AJCC pathologic stages by:
    - Converting Stage III/IIIA/IIIB/IIIC to Stage III
    - Converting Stage IV/IVA/IVB to Stage IV
    - Keeping Stage I and Stage II as is
    - Removing rows with other values
    
    Parameters:
    df (pandas.DataFrame): Input DataFrame with 'ajcc_pathologic_stage' column
    
    Returns:
    pandas.DataFrame: DataFrame with standardized stages and filtered rows
    """
    # Create a copy to avoid modifying the original DataFrame
    df_cleaned = df.copy()
    
    # Create stage mapping dictionary
    stage_mapping = {
        'Stage I': 'Stage I',
        'Stage II': 'Stage II',
        'Stage III': 'Stage III',
        'Stage IIIA': 'Stage III',
        'Stage IIIB': 'Stage III',
        'Stage IIIC': 'Stage III',
        'Stage IV': 'Stage IV',
        'Stage IVA': 'Stage IV',
        'Stage IVB': 'Stage IV'
    }
    
    # Apply mapping and filter rows
    df_cleaned['ajcc_pathologic_stage'] = df_cleaned['ajcc_pathologic_stage'].map(stage_mapping)
    df_cleaned = df_cleaned[df_cleaned['ajcc_pathologic_stage'].notna()]
    
    return df_cleaned

# Read the input file
df = pd.read_csv('clinical_fine.tsv', sep='\t')

# Apply standardization
df_standardized = standardize_pathologic_stages(df)

# Save to output file
df_standardized.to_csv('clinical_stage.tsv', sep='\t', index=False)

# Print summary of the transformation
print("\nOriginal number of rows:", len(df))
print("Number of rows after standardization:", len(df_standardized))
print("\nStandardized stage counts:")
print(df_standardized['ajcc_pathologic_stage'].value_counts())