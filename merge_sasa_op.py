import csv

# File paths
file1_path = "/dors/wankowicz_lab/ellas/apo_op/close_resi_filtered/binding_combined.csv"
file2_path = "/dors/wankowicz_lab/ellas/apo_sasa/combined_solvent_data.csv"
output_path = "/dors/wankowicz_lab/ellas/apo_sasa/sasa_op_merged_binding.csv"

# Load SASA data using (chain, resi) as key
file2_data = {}
with open(file2_path, mode='r') as file2:
    reader = csv.DictReader(file2)
    for row in reader:
        key = (row['Chain'], row['resi'])  # Ensure these column names match the SASA file
        file2_data[key] = row  # Keep all columns including pdb_id, SASA, solvent_exposure, etc.

# Merge binding OP data with SASA data
with open(file1_path, mode='r') as file1, open(output_path, mode='w', newline='') as output_file:
    reader = csv.DictReader(file1)

    # Use the first row from file2_data to determine new SASA columns
    sample_file2_row = next(iter(file2_data.values()))
    new_columns = [col for col in sample_file2_row.keys() if col not in reader.fieldnames]

    # Final output headers
    combined_headers = reader.fieldnames + new_columns
    writer = csv.DictWriter(output_file, fieldnames=combined_headers)
    writer.writeheader()

    for row in reader:
        key = (row['chain'], row['resi'])
        if key in file2_data:
            merged_row = {**row, **file2_data[key]}
        else:
            # Fill with blanks if no match
            merged_row = {**row, **{col: '' for col in new_columns}}
        writer.writerow(merged_row)

print(f"Merged data saved to {output_path}")

# Show first few rows for verification
with open(output_path, mode='r') as result_file:
    result_reader = csv.DictReader(result_file)
    print("\nFirst few rows of the resulting CSV:")
    for i, line in enumerate(result_reader):
        print(line)
        if i >= 4:
            break
