import argparse
import csv
import os
import glob


def gather_rfree_values(input_csv, base_path, intermediate_csv, final_csv, log_file):
    rfree_data = {}
    log_messages = []
    single_pdb_map = {}

    # Step 1: Preload all file paths for faster searching
    print("Preloading all file paths...", flush=True)
    folders = [
        "10001_20000", "20001_30000", "40001_50000", "60001_70000", "80001_end",
        "1_10000", "30001_40000", "50001_60000", "70001_80000"
    ]
    all_files = []
    for folder in folders:
        search_path = os.path.join(base_path, folder, "output", "*_rvalues.csv")
        all_files.extend(glob.glob(search_path))
    file_map = {os.path.basename(f).split('_rvalues.csv')[0]: f for f in all_files}
    print(f"Total files preloaded: {len(file_map)}", flush=True)

    # Step 2: Read the input CSV
    log_messages.append(f"Attempting to read input CSV: {input_csv}")
    try:
        with open(input_csv, mode='r') as infile:
            reader = csv.DictReader(infile)
            uniprot_pdb_map = {}
            row_count = 0
            for row in reader:
                row_count += 1
                if row_count % 1000 == 0:
                    print(f"Processed {row_count} rows...", flush=True)
                try:
                    uniprot_id = row['uniprot_id']
                    pdb_ids = row['pdb_ids'].split(',')
                    if uniprot_id not in uniprot_pdb_map:
                        uniprot_pdb_map[uniprot_id] = set()
                    for pdb_id in pdb_ids:
                        uniprot_pdb_map[uniprot_id].add(pdb_id.strip())
                except KeyError:
                    log_messages.append("Error: Unable to extract 'uniprot_id' and 'pdb_ids' from the row")
                    continue
            log_messages.append(f"Total rows processed: {row_count}")
    except Exception as e:
        log_messages.append(f"Error reading input CSV: {str(e)}")
        with open(log_file, mode='w') as logfile:
            for message in log_messages:
                logfile.write(message + '\n')
        return

    # Step 3: Gather Rfree values
    print("Gathering Rfree values...", flush=True)
    for uniprot_id, pdb_ids in uniprot_pdb_map.items():
        pdb_ids = list(pdb_ids)
        if len(pdb_ids) == 1:
            single_pdb_map[uniprot_id] = pdb_ids[0]
            continue
        rfree_data[uniprot_id] = []
        for pdb_id in pdb_ids:
            if pdb_id not in file_map:
                log_messages.append(f"File not found for PDB ID {pdb_id} under UniProt ID {uniprot_id}")
                continue

            file_path = file_map[pdb_id]
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

    # Step 4: Write intermediate CSV
    print("Writing intermediate CSV...", flush=True)
    with open(intermediate_csv, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['UniProtID', 'PDBID', 'Rfree'])
        for uniprot_id, values in rfree_data.items():
            for pdb_id, rfree_value in values:
                writer.writerow([uniprot_id, pdb_id, rfree_value])

    # Step 5: Find the best PDB ID per UniProt ID
    print("Finding the best PDB IDs...", flush=True)
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
    print("Script completed successfully.", flush=True)

    # Step 6: Log missing files and errors
    with open(log_file, mode='w') as logfile:
        for message in log_messages:
            logfile.write(message + '\n')


if __name__ == '__main__':
    input_csv = "/dors/wankowicz_lab/ellas/holo_uniprot_summary.csv"
    base_path = "/dors/wankowicz_lab/all_pdb/"
    intermediate_csv = "/dors/wankowicz_lab/ellas/holo_uniprot_summary_intermediate.csv"
    final_csv = "/dors/wankowicz_lab/ellas/holo_uniprot_summary_bestrfree_pdbs.csv"
    log_file = "/dors/wankowicz_lab/ellas/holo_uniprot_summary_rfree_log.txt"
    gather_rfree_values(input_csv, base_path, intermediate_csv, final_csv, log_file)
