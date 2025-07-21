import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.api.types import CategoricalDtype
from scipy.stats import kruskal
import logging

# === SETUP LOGGING ===
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# === INPUTS ===
pdb_list_path = "/.../updated_main_dataset.txt"
pdb_dir = "/dors/wankowicz_lab/all_pdb"
binding_residue_csv = "/.../updated_all_binding_residues.csv"
rigidity_csv = "/.../updated_avg_binding_site_op_per_pdb_with_rigidity_category.csv"
output_dir = "/.../water_analysis_outputs"
os.makedirs(output_dir, exist_ok=True)

subdirs = [
    "1_10000", "10001_20000", "20001_30000", "30001_40000",
    "40001_50000", "50001_60000", "60001_70000", "70001_80000", "80001_end"
]

# === LOAD CSVs ===
with open(pdb_list_path) as f:
    pdb_ids = [line.strip().lower() for line in f if line.strip()]

rigid_df = pd.read_csv(rigidity_csv)
rigid_df["pdb_id"] = rigid_df["pdb_id"].str.lower().str.strip()

binding_df = pd.read_csv(binding_residue_csv)
binding_df["pdb_id"] = binding_df["pdb_id"].astype(str).str.lower().str.strip()
binding_df["chain"] = binding_df["chain"].astype(str).str.strip()
binding_df["residue_label"] = binding_df["residue_label"].astype(str).str.strip()

# === OUTPUT CONTAINERS ===
all_waters = []
waters_per_pdb = []
binding_water_summary = []
found_pdbs = 0

parser = PDBParser(QUIET=True)

for i, pdb_id in enumerate(pdb_ids):
    if i % 500 == 0:
        logging.info(f"Processing PDB {i}/{len(pdb_ids)}: {pdb_id}")
    pdb_path = None
    upper_id = pdb_id.upper()

    for sub in subdirs:
        pdb_folder = os.path.join(pdb_dir, sub, upper_id)
        f1 = os.path.join(pdb_folder, f"{upper_id}.pdb")
        f2 = os.path.join(pdb_folder, f"{upper_id}_updated.pdb")
        if os.path.isfile(f1):
            pdb_path = f1
            break
        elif os.path.isfile(f2):
            pdb_path = f2
            break

    if pdb_path is None:
        logging.warning(f"PDB file not found for {pdb_id}")
        continue

    found_pdbs += 1

    try:
        structure = parser.get_structure(pdb_id, pdb_path)
        model = structure[0]

        # Count all residues
        residues = {(res.get_parent().id, res.id[1]) for res in model.get_residues() if res.id[0] == ' '}
        num_total_residues = len(residues)

        # Count water atoms
        water_atoms = [atom for atom in model.get_atoms() if atom.get_parent().get_resname() in ("HOH", "WAT")]
        waters_per_pdb.append({
            "pdb_id": pdb_id,
            "num_waters": len(water_atoms),
            "num_residues": num_total_residues
        })

        for atom in water_atoms:
            res = atom.get_parent()
            all_waters.append({
                "pdb_id": pdb_id,
                "chain": res.get_parent().id,
                "resname": res.get_resname(),
                "resid": res.id[1],
                "x": atom.coord[0],
                "y": atom.coord[1],
                "z": atom.coord[2]
            })

        # Nearby water analysis for binding residues
        binding_res = binding_df[binding_df["pdb_id"] == pdb_id]
        if binding_res.empty or len(water_atoms) == 0:
            continue

        chain_res_map = {(row["chain"], row["residue_label"]) for _, row in binding_res.iterrows()}
        binding_atoms = [
            atom for chain in model
            for res in chain
            for atom in res
            if (chain.id, str(res.id[1])) in chain_res_map
        ]

        if not binding_atoms:
            continue

        ns = NeighborSearch(water_atoms)
        all_close = []
        for atom in binding_atoms:
            close = ns.search(atom.coord, 5.0, level="A")
            all_close.extend(close)

        total_close = len(all_close)
        unique_close = len(set((a.get_parent().id, a.get_parent().get_parent().id) for a in all_close))

        binding_water_summary.append({
            "pdb_id": pdb_id,
            "num_binding_residues": len(chain_res_map),
            "total_nearby_waters": total_close,
            "unique_nearby_waters": unique_close
        })

    except Exception as e:
        logging.error(f"Error processing {pdb_id}: {e}")

logging.info(f"Found valid PDBs for {found_pdbs} / {len(pdb_ids)} structures.")

# === SAVE CSVs ===
pd.DataFrame(all_waters).to_csv(os.path.join(output_dir, "updated_all_water_atoms.csv"), index=False)
df_total = pd.DataFrame(waters_per_pdb)
df_total.to_csv(os.path.join(output_dir, "updated_water_count_per_pdb.csv"), index=False)
df_binding = pd.DataFrame(binding_water_summary)
df_binding.to_csv(os.path.join(output_dir, "updated_binding_site_water_summary.csv"), index=False)

# === MERGE ===
merged_total = pd.merge(df_total, rigid_df, on="pdb_id")
merged_binding = pd.merge(df_binding, rigid_df, on="pdb_id")

# === ADD NORMALIZED COLUMNS ===
merged_total["waters_per_residue"] = merged_total["num_waters"] / merged_total["num_residues"]
merged_binding["waters_per_binding_residue"] = merged_binding["unique_nearby_waters"] / merged_binding["num_binding_residues"]

# === RIGIDITY ORDER ===
cat_order = ["very flexible", "flexible", "rigid", "very rigid"]
rigidity_dtype = CategoricalDtype(categories=cat_order, ordered=True)
for df in [merged_total, merged_binding]:
    df["rigidity_quartile"] = df["rigidity_quartile"].astype(rigidity_dtype)

# === REVERSED VIRIDIS-LIKE PALETTE (very rigid = darkest purple) ===
color_map = {
    "very rigid": "#440154",     # darkest
    "rigid": "#31688e",
    "flexible": "#35b779",
    "very flexible": "#fde725"   # brightest
}

def plot_box(df, ycol, filename, ylabel, remove_outliers=False):
    data = df.copy()
    if remove_outliers:
        def rm_outliers(grp):
            q1 = grp[ycol].quantile(0.25)
            q3 = grp[ycol].quantile(0.75)
            iqr = q3 - q1
            return grp[(grp[ycol] >= q1 - 1.5 * iqr) & (grp[ycol] <= q3 + 1.5 * iqr)]
        data = data.groupby("rigidity_quartile", group_keys=False).apply(rm_outliers)

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=data, x="rigidity_quartile", y=ycol, palette=color_map, order=cat_order)
    plt.xlabel("")  # remove x-axis label
    plt.xticks(ticks=range(4), labels=[label.title() for label in cat_order], fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel(ylabel, fontsize=18)
    plt.tight_layout()
    suffix = "_no_outliers" if remove_outliers else ""
    plt.savefig(os.path.join(output_dir, f"{filename}{suffix}.png"), dpi=300)
    plt.close()

    # Kruskal–Wallis test
    groups = [grp[ycol].dropna().values for _, grp in data.groupby("rigidity_quartile")]
    if len(groups) >= 2:
        stat, p = kruskal(*groups)
        logging.info(f"Kruskal–Wallis for {filename}: H={stat:.3f}, p={p:.3e}")
    else:
        logging.info(f"Insufficient groups for Kruskal–Wallis: {filename}")

# === PLOT ALL ===
plot_box(merged_total, "num_waters", "boxplot_total_waters_per_pdb", "Total Waters per PDB")
plot_box(merged_binding, "unique_nearby_waters", "boxplot_binding_site_nearby_waters", "Unique Waters Near Binding Site (≤5 Å)")
plot_box(merged_total, "waters_per_residue", "boxplot_total_normalized", "Waters per Residue")
plot_box(merged_binding, "waters_per_binding_residue", "boxplot_binding_normalized", "Waters per Binding Site Residue")

# === OUTLIER-REMOVED ===
plot_box(merged_total, "num_waters", "boxplot_total_waters_per_pdb", "Total Waters per PDB", remove_outliers=True)
plot_box(merged_binding, "unique_nearby_waters", "boxplot_binding_site_nearby_waters", "Unique Waters Near Binding Site (≤5 Å)", remove_outliers=True)
plot_box(merged_total, "waters_per_residue", "boxplot_total_normalized", "Waters per Residue", remove_outliers=True)
plot_box(merged_binding, "waters_per_binding_residue", "boxplot_binding_normalized", "Waters per Binding Site Residue", remove_outliers=True)

logging.info("Water analysis script finished.")
