import csv
import shutil
import os

def copy_sra_files(csv_file, source_dir, dest_dir):
    # Ensure the destination directory exists
    os.makedirs(dest_dir, exist_ok=True)

    # Read the CSV file and extract SRA IDs
    sra_ids = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            sra_ids.append(row['sra_id'])

    # Copy files
    for sra_id in sra_ids:
        source_file = os.path.join(source_dir, sra_id, f"{sra_id}.tsv")
        dest_file = os.path.join(dest_dir, f"{sra_id}.tsv")
        
        try:
            shutil.copy2(source_file, dest_file)
            print(f"Copied {sra_id}.tsv successfully")
        except FileNotFoundError:
            print(f"Error: File {source_file} not found")
        except PermissionError:
            print(f"Error: Permission denied when copying {sra_id}.tsv")
        except Exception as e:
            print(f"Error copying {sra_id}.tsv: {str(e)}")

# Main execution
csv_file = 'all_family.csv'
source_dir = '/mnt/d'
dest_dir = '.'  # Current directory

copy_sra_files(csv_file, source_dir, dest_dir)

print("File copying process completed.")