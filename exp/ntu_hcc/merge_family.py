import csv

def read_csv(filename):
    data = {}
    with open(filename, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            data[row['sample_id']] = row
    return data

def read_id_list(filename):
    id_map = {}
    with open(filename, 'r') as file:
        for line in file:
            gsm_id, sra_id = line.strip().split()
            id_map[gsm_id] = sra_id
    return id_map

def merge_data(hcc_data, liv_data, hcc_id_map, liv_id_map):
    merged_data = []
    
    # Process HCC data
    for gsm_id, row in hcc_data.items():
        merged_row = {
            'sample': hcc_id_map.get(gsm_id, ''),
            'etiology': row['etiology'],
            'ctnnb1 status': row['ctnnb1 status'],
            'efs days': row['efs days'],
            'os days': row['os days']
        }
        merged_data.append(merged_row)
    
    # Process Liver data
    for gsm_id, row in liv_data.items():
        merged_row = {
            'sample': liv_id_map.get(gsm_id, ''),
            'etiology': row['disease'],
            'ctnnb1 status': row['disease'],
            'efs days': '',
            'os days': ''
        }
        merged_data.append(merged_row)
    
    return merged_data

def write_csv(data, filename):
    fieldnames = ['sample', 'etiology', 'ctnnb1 status', 'efs days', 'os days']
    with open(filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

# Main execution
hcc_file = 'hcc_family.csv'
liv_file = 'liv_family.csv'
hcc_id_file = 'hcc_list.txt'
liv_id_file = 'liv_list.txt'
output_file = 'all_family.csv'

hcc_data = read_csv(hcc_file)
liv_data = read_csv(liv_file)
hcc_id_map = read_id_list(hcc_id_file)
liv_id_map = read_id_list(liv_id_file)

merged_data = merge_data(hcc_data, liv_data, hcc_id_map, liv_id_map)

write_csv(merged_data, output_file)

print(f"Merged data has been written to {output_file}")