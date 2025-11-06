#!/usr/bin/env python3
"""
Parse all FreeSASA JSON outputs into CSV format with per-residue statistics.
Usage: python parse_freesasa_json.py
"""

import json
import csv
import os
from glob import glob
from datetime import datetime

# === CONFIG ===
JSON_DIR = "/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_per_atom"
CSV_DIR = "/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_per_atom_csv"
LOG_FILE = "/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/parse_freesasa_updated.log"

def log_message(message):
    """Log message with timestamp"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)
    with open(LOG_FILE, 'a') as f:
        f.write(log_msg + '\n')

def parse_freesasa_json(json_file, output_csv):
    """
    Parse FreeSASA JSON output and create CSV with per-residue data.
    
    Columns:
    - pdb_id
    - chain
    - resi (residue number)
    - resn (residue name)
    - total_sasa
    - apolar_sasa
    - polar_sasa
    - n_atoms_total
    - n_atoms_apolar
    - n_atoms_polar
    - total_sasa_normalized (by total atoms)
    - apolar_sasa_normalized (by apolar atoms)
    - polar_sasa_normalized (by polar atoms)
    - is_exposed (1 if total_sasa > 0, else 0)
    """
    
    # Load JSON
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract PDB ID from filename
    pdb_id = os.path.basename(json_file).replace('_sasa_atom.json', '')
    
    # Prepare output data
    rows = []
    
    # Navigate through the JSON structure
    for result in data.get('results', []):
        for structure in result.get('structure', []):
            for chain_data in structure.get('chains', []):
                chain = chain_data.get('label', '')
                
                for residue in chain_data.get('residues', []):
                    resn = residue.get('name', '')
                    resi = residue.get('number', '')
                    
                    # Get area information
                    area = residue.get('area', {})
                    total_sasa = area.get('total', 0.0)
                    polar_sasa = area.get('polar', 0.0)
                    apolar_sasa = area.get('apolar', 0.0)
                    
                    # Count atoms by polarity
                    atoms = residue.get('atoms', [])
                    n_atoms_total = len(atoms)
                    n_atoms_polar = sum(1 for atom in atoms if atom.get('is-polar', False))
                    n_atoms_apolar = n_atoms_total - n_atoms_polar
                    
                    # Calculate normalized values
                    total_sasa_normalized = total_sasa / n_atoms_total if n_atoms_total > 0 else 0.0
                    apolar_sasa_normalized = apolar_sasa / n_atoms_apolar if n_atoms_apolar > 0 else 0.0
                    polar_sasa_normalized = polar_sasa / n_atoms_polar if n_atoms_polar > 0 else 0.0
                    
                    # Determine if exposed
                    is_exposed = 1 if total_sasa > 0 else 0
                    
                    # Append row
                    rows.append({
                        'pdb_id': pdb_id,
                        'chain': chain,
                        'resi': resi,
                        'resn': resn,
                        'total_sasa': total_sasa,
                        'apolar_sasa': apolar_sasa,
                        'polar_sasa': polar_sasa,
                        'n_atoms_total': n_atoms_total,
                        'n_atoms_apolar': n_atoms_apolar,
                        'n_atoms_polar': n_atoms_polar,
                        'total_sasa_normalized': total_sasa_normalized,
                        'apolar_sasa_normalized': apolar_sasa_normalized,
                        'polar_sasa_normalized': polar_sasa_normalized,
                        'is_exposed': is_exposed
                    })
    
    # Write to CSV
    if rows:
        fieldnames = [
            'pdb_id', 'chain', 'resi', 'resn',
            'total_sasa', 'apolar_sasa', 'polar_sasa',
            'n_atoms_total', 'n_atoms_apolar', 'n_atoms_polar',
            'total_sasa_normalized', 'apolar_sasa_normalized', 'polar_sasa_normalized',
            'is_exposed'
        ]
        
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        
        return True, len(rows)
    else:
        return False, 0

def main():
    """Process all JSON files in the input directory"""
    
    # Create output directory
    os.makedirs(CSV_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    
    log_message("Starting JSON to CSV conversion...")
    
    # Find all JSON files
    json_files = sorted(glob(os.path.join(JSON_DIR, "*_sasa_atom.json")))
    
    if not json_files:
        log_message(f"ERROR: No JSON files found in {JSON_DIR}")
        raise SystemExit(1)
    
    log_message(f"Found {len(json_files)} JSON files to process")
    
    success_count = 0
    fail_count = 0
    skip_count = 0
    
    # Process each JSON file
    for idx, json_file in enumerate(json_files, 1):
        pdb_id = os.path.basename(json_file).replace('_sasa_atom.json', '')
        output_csv = os.path.join(CSV_DIR, f"{pdb_id}_freesasa.csv")
        
        # Skip if output already exists
        if os.path.exists(output_csv):
            log_message(f"[{idx}/{len(json_files)}] ⊙ Skipping: {pdb_id} (CSV already exists)")
            skip_count += 1
            continue
        
        log_message(f"[{idx}/{len(json_files)}] Processing: {pdb_id}")
        
        try:
            success, n_residues = parse_freesasa_json(json_file, output_csv)
            if success:
                log_message(f"[{idx}/{len(json_files)}] ✔ Success: {pdb_id} ({n_residues} residues)")
                success_count += 1
            else:
                log_message(f"[{idx}/{len(json_files)}] ✗ Failed: {pdb_id} (no residues found)")
                fail_count += 1
                # Fail on first file if it produces empty output
                if idx == 1:
                    log_message("FATAL: First file failed - stopping execution")
                    raise SystemExit(1)
        except Exception as e:
            log_message(f"[{idx}/{len(json_files)}] ✗ Failed: {pdb_id} - Error: {str(e)}")
            fail_count += 1
            # Fail on first file if it errors
            if idx == 1:
                log_message("FATAL: First file failed - stopping execution")
                raise SystemExit(1)
    
    # Summary
    log_message("")
    log_message("Summary:")
    log_message(f"  - Successfully parsed: {success_count} files")
    log_message(f"  - Skipped (already exists): {skip_count} files")
    log_message(f"  - Failed: {fail_count} files")
    log_message("Complete!")

if __name__ == "__main__":
    main()
