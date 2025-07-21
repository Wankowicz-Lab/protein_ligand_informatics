import os
import pandas as pd
from Bio.PDB import PDBParser
import numpy as np

# === File paths ===
print("Loading input files...")

rigidity_file = "/.../avg_binding_site_op_per_pdb_with_rigidity_category.csv"
water_file = "/.../updated_binding_site_water_summary.csv"
hbond_file = "/.../water_bond_clash_summary_all_pdb.csv"
temp_file = "/.../pdb_temperatures.csv"
binding_residue_file = "/.../all_binding_residues.csv"

pdb_root = "/dors/wankowicz_lab/all_pdb"
subdirs = ["1_10000", "10001_20000", "20001_30000", "30001_40000", 
           "40001_50000", "50001_60000", "60001_70000", "70001_80000", "80001_end"]

# === Load metadata ===
rigidity_df = pd.read_csv(rigidity_file)
water_df = pd.read_csv(water_file)
hbond_df = pd.read_csv(hbond_file)
temp_df = pd.read_csv(temp_file)
binding_df = pd.read_csv(binding_residue_file)

print("Cleaning and merging metadata...")

# Normalize PDB IDs
# Standardize all relevant PDB ID columns
rigidity_df['pdb_id'] = rigidity_df['pdb_id'].astype(str).str.strip().str.lower()
water_df['pdb_id'] = water_df['pdb_id'].astype(str).str.strip().str.lower()
hbond_df['pdb'] = hbond_df['pdb'].astype(str).str.strip().str.lower()
temp_df['pdb_id'] = temp_df['pdb_id'].astype(str).str.strip().str.lower()
binding_df['pdb_id'] = binding_df['pdb_id'].astype(str).str.strip().str.lower()


# Merge metadata
merged = rigidity_df.merge(
    water_df[['pdb_id', 'unique_nearby_waters']], on='pdb_id', how='left'
)
merged = merged.merge(
    hbond_df[['pdb', 'total_bonds']].rename(columns={'pdb': 'pdb_id'}), on='pdb_id', how='left'
)
merged = merged.merge(
    temp_df[['pdb_id', 'temp']].rename(columns={'temp': 'temperature_K'}), on='pdb_id', how='left'
)


# Add temperature category
def classify_temp(temp):
    if pd.isna(temp):
        return "Unknown"
    elif temp >= 293:
        return "RT"
    elif temp < 200:
        return "Cryo"
    else:
        return "Other"
# Convert temperature column to numeric
merged["temperature_K"] = pd.to_numeric(merged["temperature_K"], errors="coerce")

# Assign temperature category
def classify_temp(temp):
    if pd.isna(temp):
        return "Unknown"
    elif temp >= 293:
        return "RT"
    elif temp < 200:
        return "Cryo"
    else:
        return "Other"

merged["temperature_category"] = merged["temperature_K"].apply(classify_temp)


print(f"Total PDBs for analysis: {len(merged)}")

# === Setup B-factor parser ===
parser = PDBParser(QUIET=True)
b_results = []

# Loop through PDBs
for i, row in merged.iterrows():
    pdb_id = row['pdb_id']
    pdb_id_upper = pdb_id.upper()
    found = False
    pdb_path = None

    for sub in subdirs:
        base = os.path.join(pdb_root, sub, pdb_id_upper)
        for suffix in [f"{pdb_id_upper}_updated.pdb", f"{pdb_id_upper}.pdb"]:
            path = os.path.join(base, suffix)
            if os.path.isfile(path):
                pdb_path = path
                found = True
                break
        if found:
            break

    if not found:
        print(f"[{i+1}/{len(merged)}] Missing PDB: {pdb_id}")
        b_results.append({
            "pdb_id": pdb_id,
            "avg_bfactor_binding_residues": np.nan,
            "avg_bfactor_binding_waters": np.nan
        })
        continue

    try:
        structure = parser.get_structure(pdb_id, pdb_path)
        atoms = list(structure.get_atoms())

        # Get binding site residue keys
        bind_res = binding_df[binding_df["pdb_id"] == pdb_id]
        binding_res_keys = set((row["chain"], int(row["residue_label"])) for _, row in bind_res.iterrows())

        b_res = []
        water_b = []

        for atom in atoms:
            res = atom.get_parent()
            if res.id[0] != ' ':
                continue  # Skip heteroatoms
            chain_id = res.get_parent().id
            resnum = res.id[1]
            if atom.name == "CA" and (chain_id, resnum) in binding_res_keys:
                b_res.append(atom.bfactor)

        binding_atoms = [atom for atom in atoms if
                         atom.name == "CA" and
                         (atom.get_parent().get_parent().id, atom.get_parent().id[1]) in binding_res_keys]
        water_atoms = [atom for atom in atoms if atom.get_parent().resname in ["HOH", "WAT"] and atom.name == "O"]

        for wat in water_atoms:
            for b_atom in binding_atoms:
                dist = np.linalg.norm(wat.coord - b_atom.coord)
                if dist <= 3.2:
                    water_b.append(wat.bfactor)
                    break

        b_results.append({
            "pdb_id": pdb_id,
            "avg_bfactor_binding_residues": np.mean(b_res) if b_res else np.nan,
            "avg_bfactor_binding_waters": np.mean(water_b) if water_b else np.nan
        })
        if (i + 1) % 100 == 0:
            print(f"[{i+1}/{len(merged)}] Processed {pdb_id}")

    except Exception as e:
        print(f"[{i+1}/{len(merged)}] Error parsing {pdb_id}: {e}")
        b_results.append({
            "pdb_id": pdb_id,
            "avg_bfactor_binding_residues": np.nan,
            "avg_bfactor_binding_waters": np.nan
        })

# Final merge
print("Merging results...")
b_df = pd.DataFrame(b_results)
final_df = merged.merge(b_df, on="pdb_id", how="left")

output_path = "/.../rt_cryo_rigidity_bfactor_analysis.csv"
final_df.to_csv(output_path, index=False)
print(f"Saved results to: {output_path}")
