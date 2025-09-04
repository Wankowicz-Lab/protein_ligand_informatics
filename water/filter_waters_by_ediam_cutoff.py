import os
import pandas as pd

# === INPUTS ===
pdb_list_file = "..."
output_csv = "..."
base_dir = "/dors/wankowicz_lab/PDBRedo/pdb-redo3/"

# === Read list of PDB IDs ===
with open(pdb_list_file, "r") as f:
    pdb_ids = [line.strip().lower() for line in f if len(line.strip()) == 4]

print(f"Processing {len(pdb_ids)} PDB IDs for EDIAm water filtering...")

# === Output accumulator ===
all_rows = []

for i, pdb_id in enumerate(pdb_ids):
    if i % 100 == 0:
        print(f"Progress: {i}/{len(pdb_ids)} ({i/len(pdb_ids)*100:.1f}%)")
    
    middle_two = pdb_id[1:3]
    csv_path = os.path.join(base_dir, middle_two, pdb_id, f"{pdb_id}_redo_residue_stats.csv")
    
    if not os.path.isfile(csv_path):
        print(f"[SKIP] {csv_path} not found")
        continue
    
    try:
        df = pd.read_csv(csv_path)
        
        # Check if we have the expected columns
        required_cols = ['EDIAm', 'pdb_compID', 'pdb_strandID', 'pdb_seqNum', 'pdb_insCode']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"[WARNING] {pdb_id}: Missing columns {missing_cols}")
            continue
        
        # Filter waters with EDIAm > 0.3
        water_mask = (df['pdb_compID'].isin(["HOH", "WAT"])) & (df["EDIAm"] > 0.3)
        filtered = df[water_mask].copy()
        
        if len(filtered) == 0:
            continue  # Skip logging for cleaner output
        
        # Select and rename columns
        result = filtered[required_cols].copy()
        result = result.rename(columns={
            'pdb_strandID': 'chain',
            'pdb_seqNum': 'resid',
            'pdb_insCode': 'insertion_code',
            'pdb_compID': 'resname'
        })
        
        # Add PDB ID
        result.insert(0, "pdb_id", pdb_id.upper())
        
        # Clean up insertion codes - handle NaN, None, and convert to string
        result['insertion_code'] = result['insertion_code'].fillna('').astype(str).str.strip()
        # Replace 'nan' string with empty string (sometimes happens after astype(str))
        result['insertion_code'] = result['insertion_code'].replace('nan', '')
        
        # Also ensure other columns are proper types
        result['chain'] = result['chain'].astype(str).str.strip()
        result['resid'] = result['resid'].astype(int)
        result['resname'] = result['resname'].astype(str).str.strip()
        
        all_rows.append(result)
        print(f"[DONE] {pdb_id}: {len(result)} waters with EDIAm > 0.3")
        
        # Show sample for first PDB
        if i == 0 and len(result) > 0:
            print(f"[SAMPLE] First water example: {result.iloc[0].to_dict()}")
    
    except Exception as e:
        print(f"[ERROR] Failed to process {pdb_id}: {e}")

# === Save final merged CSV ===
if all_rows:
    final_df = pd.concat(all_rows, ignore_index=True)
    final_df.to_csv(output_csv, index=False)
    
    print(f"EDIAm-filtered waters saved to: {output_csv}")
    print(f"Summary:")
    print(f"   - Total waters meeting EDIAm > 0.3: {len(final_df)}")
    print(f"   - Unique PDBs: {final_df['pdb_id'].nunique()}")
    print(f"   - EDIAm range: {final_df['EDIAm'].min():.3f} - {final_df['EDIAm'].max():.3f}")
    print(f"   - Mean EDIAm: {final_df['EDIAm'].mean():.3f}")
    # Convert chains to string for sorting
    unique_chains = [str(x) for x in final_df['chain'].unique()]
    print(f"   - Chains: {sorted(unique_chains)}")
    print(f"   - Final columns: {final_df.columns.tolist()}")
    
    # Show distribution by PDB
    pdb_counts = final_df['pdb_id'].value_counts()
    print(f"   - Waters per PDB range: {pdb_counts.min()} - {pdb_counts.max()}")
    print(f"   - Mean waters per PDB: {pdb_counts.mean():.1f}")
    
else:
    print("\n No waters found meeting EDIAm > 0.3 criteria in any PDBs.")

print("\nEDIAm water filtering complete!")
