import os
import pandas as pd

# === INPUTS ===
pdb_list_file = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset.txt"
output_csv = "/dors/wankowicz_lab/ellas/main_dataset/confident_water_EDIAm.csv"
base_dir = "/dors/wankowicz_lab/PDBRedo/pdb-redo3/"

# === Read list of PDB IDs ===
with open(pdb_list_file, "r") as f:
    pdb_ids = [line.strip().lower() for line in f if len(line.strip()) == 4]

# === Output accumulator ===
all_rows = []

for pdb_id in pdb_ids:
    middle_two = pdb_id[1:3]
    csv_path = os.path.join(base_dir, middle_two, pdb_id, f"{pdb_id}_redo_residue_stats.csv")

    if not os.path.isfile(csv_path):
        print(f"[SKIP] {csv_path} not found")
        continue

    try:
        df = pd.read_csv(csv_path)

        # Check column casing (some versions use pdb_compID or compID)
        compid_col = "pdb_compID" if "pdb_compID" in df.columns else "compID"

        # Filter waters with EDIAm > 0.3
        filtered = df[
            (df[compid_col].isin(["HOH", "WAT"])) &
            (df["EDIAm"] > 0.3)
        ][["EDIAm", compid_col]].copy()
        filtered.insert(0, "pdb_id", pdb_id)

        all_rows.append(filtered)
        print(f"[DONE] Processed {pdb_id} with {len(filtered)} confident waters")

    except Exception as e:
        print(f"[ERROR] Failed to process {pdb_id}: {e}")

# === Save final merged CSV ===
if all_rows:
    final_df = pd.concat(all_rows, ignore_index=True)
    final_df.to_csv(output_csv, index=False)
    print(f"\n utput saved to: {output_csv}")
else:
    print("\n No valid water rows found in any PDBs.")
