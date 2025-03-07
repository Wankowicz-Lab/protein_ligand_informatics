import requests
import csv
import json
import os

# Input and output file paths
input_csv = "/dors/wankowicz_lab/ellas/apo_pdb_uniprot_mapping.csv"
output_csv = "apo_merged_uniprot_data.csv"
log_file = "apo_merged_uniprot_data_log.txt"

def get_protein_family(uniprot_id):
    """Fetch protein family information for a given UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)

    if response.status_code != 200:
        log_message(f"Error: Unable to fetch data for {uniprot_id}")
        return None

    data = response.json()
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SIMILARITY":
            for entry in comment.get("texts", []):
                if "Belongs to" in entry.get("value", ""):
                    return entry["value"]
    return None

def get_enzyme_number(uniprot_id):
    """Fetch EC number for a given UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)

    if response.status_code != 200:
        log_message(f"Error: Unable to fetch data for {uniprot_id}")
        return None

    data = response.json()

    # Check if EC number exists under recommendedName
    ec_numbers = []
    protein_names = data.get("proteinDescription", {}).get("recommendedName", {})
    if "ecNumbers" in protein_names:
        ec_numbers = [ec["value"] for ec in protein_names["ecNumbers"]]

    # If no EC number is found, check the enzyme and pathway databases
    if not ec_numbers:
        for db_entry in data.get("enzymeAndPathway", []):
            if db_entry.get("database") == "BRENDA" or db_entry.get("database") == "ENZYME":
                ec_numbers.append(db_entry.get("id"))

    return "; ".join(ec_numbers) if ec_numbers else None

def log_message(message):
    """Log messages to a file."""
    with open(log_file, "a") as log:
        log.write(message + "\n")
    print(message)

# Process the CSV file
with open(input_csv, "r") as infile, open(output_csv, "w", newline="") as outfile:
    reader = csv.DictReader(infile)
    fieldnames = ["pdb_id", "uniprot_id", "protein_family", "enzyme_number"]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        pdb_id = row["pdb_id"]
        uniprot_id = row["uniprot_id"]

        log_message(f"Processing UniProt ID: {uniprot_id}")

        protein_family = get_protein_family(uniprot_id)
        enzyme_number = get_enzyme_number(uniprot_id)

        if not protein_family and not enzyme_number:
            log_message(f"Skipping {uniprot_id}: No protein family or enzyme number found.")

        writer.writerow({
            "pdb_id": pdb_id,
            "uniprot_id": uniprot_id,
            "protein_family": protein_family if protein_family else "",
            "enzyme_number": enzyme_number if enzyme_number else ""
        })

log_message("Processing complete. Output saved to output_uniprot_data.csv")
