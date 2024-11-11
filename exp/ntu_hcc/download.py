import requests
import json
import os
from urllib.request import urlretrieve

def explore_available_files():
    """Explore what RNA-seq files are available for TCGA-LIHC"""
    
    files_endpt = "https://api.gdc.cancer.gov/files"
    
    # Set up the query filters - just looking for TCGA-LIHC RNA-seq files
    filters = {
        "op": "and",
        "content":[
            {
                "op": "=",
                "content":{
                    "field": "cases.project.project_id",
                    "value": "TCGA-LIHC"
                }
            },
            {
                "op": "=",
                "content":{
                    "field": "data_type",
                    "value": "Gene Expression Quantification"
                }
            }
        ]
    }

    # Define the query parameters
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,data_type,data_format,experimental_strategy",
        "format": "JSON",
        "size": "1000"
    }
    
    # Make the API call
    response = requests.get(files_endpt, params=params)
    
    # Parse the results
    file_data = json.loads(response.content.decode("utf-8"))
    
    return file_data["data"]["hits"]

def get_file_ids(file_name_pattern):
    """Get file IDs for RNA-seq gene counts from TCGA-LIHC"""
    
    files_endpt = "https://api.gdc.cancer.gov/files"
    
    # Set up the query filters
    filters = {
        "op": "and",
        "content":[
            {
                "op": "=",
                "content":{
                    "field": "cases.project.project_id",
                    "value": "TCGA-LIHC"
                }
            },
            {
                "op": "=",
                "content":{
                    "field": "data_type",
                    "value": "Gene Expression Quantification"
                }
            },
            {
                "op": "like",
                "content":{
                    "field": "file_name",
                    "value": f"*{file_name_pattern}*"
                }
            }
        ]
    }

    # Define the query parameters
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.submitter_id",
        "format": "JSON",
        "size": "1000"
    }
    
    # Make the API call
    response = requests.get(files_endpt, params=params)
    
    # Parse the results
    file_data = json.loads(response.content.decode("utf-8"))
    
    return file_data["data"]["hits"]

def download_files(file_hits, output_dir="LIHC_RNA_counts"):
    """Download files from GDC using file IDs"""
    
    data_endpt = "https://api.gdc.cancer.gov/data/"
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Download each file
    for hit in file_hits:
        file_id = hit["file_id"]
        case_id = hit["cases"][0]["submitter_id"]
        file_name = hit["file_name"]
        
        # Create the output filename
        output_file = os.path.join(output_dir, f"{case_id}_{file_name}")
        
        # Download the file
        print(f"Downloading {case_id}...")
        response = requests.get(f"{data_endpt}{file_id}", headers={"Content-Type": "application/json"})
        
        # Save the file
        with open(output_file, 'wb') as f:
            f.write(response.content)
            
        print(f"Saved to {output_file}")

def main():
    # First, let's explore what files are available
    print("Exploring available RNA-seq files in TCGA-LIHC...")
    available_files = explore_available_files()
    
    print(f"\nFound {len(available_files)} RNA-seq files. Here are some examples:")
    # Print first few unique file names to see the pattern
    unique_names = set()
    for hit in available_files[:10]:
        unique_names.add(hit["file_name"])
    
    print("\nAvailable file patterns:")
    for name in unique_names:
        print(f"- {name}")
        
    # Now try to get the specific files
    file_name_pattern = "star_gene_counts"  # Modified pattern
    print(f"\nSearching for files matching pattern: {file_name_pattern}")
    file_hits = get_file_ids(file_name_pattern)
    
    if not file_hits:
        print("No files found matching the exact pattern.")
        return
        
    print(f"\nFound {len(file_hits)} matching files")
    
    # Ask user before downloading
    response = input("\nWould you like to download these files? (yes/no): ")
    if response.lower() == 'yes':
        print("\nStarting downloads...")
        download_files(file_hits)
        print("Download complete!")
    else:
        print("Download cancelled.")

if __name__ == "__main__":
    main()