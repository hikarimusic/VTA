import os
import shutil
import pandas as pd
from collections import defaultdict

def read_n_ambiguous(filepath):
    """Read N_ambiguous value from TCGA RNA count file."""
    try:
        # Read only the first few lines to find N_ambiguous
        with open(filepath) as f:
            for line in f:
                if 'N_ambiguous' in line:
                    # Split the line and get the unstranded value (first numeric column)
                    return int(line.split('\t')[3])
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return 0

def get_case_id(filename):
    """Extract and format case ID from filename."""
    # Extract case ID (e.g., TCGA-DD-A3A6 from filename)
    case_id = '-'.join(filename.split('_')[0].split('-')[:3])
    # Replace hyphens with underscores
    return case_id.replace('-', '_')

def process_tcga_files(input_dir, output_dir):
    """Process TCGA RNA count files and copy selected files to output directory."""
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Group files by case ID
    case_files = defaultdict(list)
    for filename in os.listdir(input_dir):
        if filename.endswith('.tsv'):
            case_id = '_'.join(filename.split('_')[0].split('-')[:3])
            case_files[case_id].append(filename)
    
    # Process each case
    for case_id, files in case_files.items():
        # Create new filename
        new_filename = f"{case_id}.tsv"
        
        if len(files) == 1:
            # If only one file exists, copy and rename it
            src = os.path.join(input_dir, files[0])
            dst = os.path.join(output_dir, new_filename)
            shutil.copy2(src, dst)
            print(f"Copied and renamed file for {case_id}")
        else:
            # If multiple files exist, compare N_ambiguous values
            file_values = []
            for filename in files:
                filepath = os.path.join(input_dir, filename)
                n_ambiguous = read_n_ambiguous(filepath)
                file_values.append((filename, n_ambiguous))
            
            # Sort by N_ambiguous value and get the file with highest value
            file_values.sort(key=lambda x: x[1], reverse=True)
            selected_file = file_values[0][0]
            
            # Copy and rename the selected file
            src = os.path.join(input_dir, selected_file)
            dst = os.path.join(output_dir, new_filename)
            shutil.copy2(src, dst)
            print(f"Selected file for {case_id} with N_ambiguous = {file_values[0][1]}")

# Usage
input_directory = "LIHC_RNA_counts"
output_directory = "LIHC_RNA_counts_selected"

process_tcga_files(input_directory, output_directory)