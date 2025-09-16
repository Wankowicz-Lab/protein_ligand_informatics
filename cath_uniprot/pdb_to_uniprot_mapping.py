import requests
import pandas as pd
import time

# === INPUT FILE ===
pdb_list_file = "/.../pdb_list.txt"   # one PDB ID per line
output_csv = "/.../pdb_to_uniprot_mappings.csv"


# === LOAD PDB IDs ===
with open(pdb_list_file, "r") as f:
    pdb_ids = [line.strip().upper() for line in f if line.strip()]  # PDB IDs are typically uppercase

print(f"Loaded {len(pdb_ids)} PDB IDs")

# === QUERY UNIPROT API ===
base_url = "https://rest.uniprot.org/uniprotkb/search"
results = []

for i, pdb_id in enumerate(pdb_ids):
    print(f"Processing {i+1}/{len(pdb_ids)}: {pdb_id}")
    
    # First, let's verify this PDB ID exists by checking RCSB PDB
    try:
        pdb_check_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        pdb_response = requests.get(pdb_check_url, timeout=10)
        if pdb_response.status_code == 404:
            print(f"  PDB ID {pdb_id} does not exist in RCSB PDB database")
            results.append({
                "pdb_id": pdb_id,
                "uniprot_accession": "PDB_NOT_FOUND",
                "uniprot_id": None,
                "protein_name": None,
                "organism": None
            })
            continue
    except Exception as e:
        print(f"  Could not verify PDB ID {pdb_id}: {e}")
    
    # Try multiple query formats for PDB cross-references
    query_formats = [
        f"xref:pdb-{pdb_id}",           # Cross-reference format
        f"database:PDB {pdb_id}",       # Database format with space
        f"database:pdb {pdb_id}",       # Database format lowercase
        f"pdb:{pdb_id}",               # Simple PDB format
    ]
    
    found_results = False
    
    # Try different query formats
    for query_format in query_formats:
        params = {
            "format": "tsv",
            "query": query_format,
            "fields": "accession,id,protein_name,organism_name",
            "size": "500"  # Limit results per PDB ID
        }
        
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            
            lines = response.text.strip().split("\n")
            
            # Check if we have results (more than just header)
            if len(lines) > 1 and any(line.strip() for line in lines[1:]):
                print(f"  Found {len([l for l in lines[1:] if l.strip()])} entries with query: {query_format}")
                found_results = True
                # Skip header and process data lines
                for line in lines[1:]:
                    if line.strip():  # Skip empty lines
                        parts = line.split("\t")
                        if len(parts) >= 4:
                            accession, uniprot_id, protein_name, organism = parts[:4]
                            results.append({
                                "pdb_id": pdb_id,
                                "uniprot_accession": accession,
                                "uniprot_id": uniprot_id,
                                "protein_name": protein_name,
                                "organism": organism
                            })
                break  # Stop trying other formats if we found results
                
        except requests.exceptions.HTTPError as e:
            if response.status_code == 400:
                continue  # Try next query format
            else:
                print(f"HTTP Error with query '{query_format}': {e}")
                break
        except Exception as e:
            print(f"Error with query '{query_format}': {e}")
            continue
    
    # If no results found with any query format
    if not found_results:
        print(f"  No UniProt entries found for {pdb_id} (tried {len(query_formats)} query formats)")
        results.append({
            "pdb_id": pdb_id,
            "uniprot_accession": None,
            "uniprot_id": None,
            "protein_name": None,
            "organism": None
        })
    
    # Add a small delay to be respectful to the API
    time.sleep(0.1)

# === SAVE OUTPUT ===
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False)
print(f"\nSaved {len(results)} results to {output_csv}")
print(f"Found UniProt mappings for {len([r for r in results if r['uniprot_accession'] is not None and r['uniprot_accession'] != 'PDB_NOT_FOUND'])} PDB IDs")

