import argparse
import csv
import os
import glob


def gather_rfree_values(input_csv, base_path, intermediate_csv, final_csv, log_file):
    rfree_data = {}
    log_messages = []
    single_pdb_map = {}

    # Step 1: Read the input CSV
    log_messages.append(f"Attempting to read input CSV: {input_csv}")
    try:
        with open(input_csv, mode='r') as infile:
            reader = csv.DictReader(infile)
            uniprot_pdb_map = {}
            # Debug: Print fieldnames
            log_messages.append(f"Fieldnames found: {reader.fieldnames}")
            # Extract only necessary columns: 'uniprot_id' and 'pdb_ids'
            row_count = 0
            for row in reader:
                row_count += 1
                try:
                    uniprot_id = row['uniprot_id']
                    pdb_ids = row['pdb_ids'].split(',')  # Split comma-separated PDB IDs
                    if uniprot_id not in uniprot_pdb_map:
                        uniprot_pdb_map[uniprot_id] = set()
                    for pdb_id in pdb_ids:
                        uniprot_pdb_map[uniprot_id].add(pdb_id.strip())  # Remove any extra spaces
                except KeyError:
                    log_messages.append("Error: Unable to extract 'uniprot_id' and 'pdb_ids' from the row")
                    continue
            log_messages.append(f"Total rows processed: {row_count}")
            log_messages.append(f"Total UniProt IDs found: {len(uniprot_pdb_map)}")
    except Exception as e:
        log_messages.append(f"Error reading input CSV: {str(e)}")
        with open(log_file, mode='w') as logfile:
            for message in log_messages:
                logfile.write(message + '\n')
        return

    # Step 2: Gather Rfree values
    for uniprot_id, pdb_ids in uniprot_pdb_map.items():
        pdb_ids = list(pdb_ids)  # Convert set to list
        if len(pdb_ids) == 1:
            single_pdb_map[uniprot_id] = pdb_ids[0]
            continue
        rfree_data[uniprot_id] = []
        for pdb_id in pdb_ids:
            pattern = os.path.join(base_path, "**", "output", f"{pdb_id}_rvalues.csv")
            matching_files = glob.glob(pattern, recursive=True)

            if not matching_files:
                log_messages.append(f"File not found for PDB ID {pdb_id} under UniProt ID {uniprot_id}")
                continue

            for file_path in matching_files:
                try:
                    with open(file_path, mode='r') as csvfile:
                        csv_reader = csv.DictReader(csvfile)
                        if 'Rfree' not in csv_reader.fieldnames:
                            log_messages.append(f"Rfree column missing in {file_path}")
                            continue
                        for row in csv_reader:
                            if 'Rfree' in row:
                                rfree_value = float(row['Rfree'])
                                rfree_data[uniprot_id].append((pdb_id, rfree_value))
                except Exception as e:
                    log_messages.append(f"Error processing file {file_path}: {str(e)}")

    # Step 3: Write intermediate CSV
    with open(intermediate_csv, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['UniProtID', 'PDBID', 'Rfree'])
        for uniprot_id, values in rfree_data.items():
            for pdb_id, rfree_value in values:
                writer.writerow([uniprot_id, pdb_id, rfree_value])

    # Step 4: Find the best PDB ID per UniProt ID
    with open(final_csv, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['UniProtID', 'BestPDBID(s)', 'BestRfree'])
        for uniprot_id, values in rfree_data.items():
            if uniprot_id in single_pdb_map:
                writer.writerow([uniprot_id, single_pdb_map[uniprot_id], ''])
                continue
            if not values:
                writer.writerow([uniprot_id, '', ''])
                continue
            max_rfree = max(values, key=lambda x: x[1])[1]
            best_pdbs = [pdb_id for pdb_id, rfree in values if rfree == max_rfree]
            writer.writerow([uniprot_id, ', '.join(best_pdbs), max_rfree])

    # Step 5: Log missing files and errors
    with open(log_file, mode='w') as logfile:
        for message in log_messages:
            logfile.write(message + '\n')


if __name__ == '__main__':
    input_csv = "/dors/wankowicz_lab/ellas/uniprot_summary.csv"
    base_path = "/dors/wankowicz_lab/all_pdb/"
    intermediate_csv = "/dors/wankowicz_lab/ellas/apo_uniprot_Intermediate.csv"
    final_csv = "/dors/wankowicz_lab/ellas/apo_uniprot_bestrfree_pdbs.csv"
    log_file = "/dors/wankowicz_lab/ellas/apo_uniprot_rfree_log.txt"
    gather_rfree_values(input_csv, base_path, intermediate_csv, final_csv, log_file)
