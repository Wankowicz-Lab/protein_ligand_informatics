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
pdb_list_path = "/.../pdb_list.txt"
pdb_dir = "/dors/wankowicz_lab/PDBRedo/pdb-redo3"  # Updated base directory
binding_residue_csv = "/.../binding_site_residues.csv"
rigidity_csv = "/.../binding_site_pocket_classified_by_rigidity_quartile.csv"
ediam_water_csv = "/.../water_meeting_ediam_cutoffs.csv"
output_dir = "/.../water_analysis_outputs"
os.makedirs(output_dir, exist_ok=True)

# === LOAD CSVs ===
with open(pdb_list_path) as f:
    pdb_ids = [line.strip().lower() for line in f if line.strip()]

rigid_df = pd.read_csv(rigidity_csv)
rigid_df["pdb_id"] = rigid_df["pdb_id"].str.lower().str.strip()

binding_df = pd.read_csv(binding_residue_csv)
binding_df["pdb_id"] = binding_df["pdb_id"].astype(str).str.lower().str.strip()
binding_df["chain"] = binding_df["chain"].astype(str).str.strip()
binding_df["resi"] = binding_df["resi"].astype(str).str.strip()

# Load EDIAm-filtered waters
logging.info("Loading EDIAm-filtered water data...")
ediam_df = pd.read_csv(ediam_water_csv)
ediam_df["pdb_id"] = ediam_df["pdb_id"].astype(str).str.lower().str.strip()
ediam_df["chain"] = ediam_df["chain"].astype(str).str.strip()
# Handle empty insertion_code column - fill with empty string if NaN
ediam_df["insertion_code"] = ediam_df["insertion_code"].fillna("").astype(str).str.strip()
logging.info(f"Loaded {len(ediam_df)} EDIAm-filtered waters")
logging.info(f"EDIAm water columns: {ediam_df.columns.tolist()}")

# === OUTPUT CONTAINERS ===
all_waters_detailed = []  # Detailed water information
pocket_water_metrics = []  # Per pocket (chain+pdb) metrics
found_pdbs = 0

parser = PDBParser(QUIET=True)

def get_pdb_path(pdb_id):
    """
    Construct PDB file path using the new directory structure:
    /dors/wankowicz_lab/PDBRedo/pdb-redo3/$middletwodigitsofpdbid/$pdbid/$pdbid_final.pdb
    """
    pdb_id_lower = pdb_id.lower()
    
    # Extract middle two digits (positions 1 and 2, 0-indexed)
    if len(pdb_id_lower) >= 4:
        middle_two = pdb_id_lower[1:3]
    else:
        # Handle edge case for shorter PDB IDs (shouldn't happen with standard PDB format)
        logging.warning(f"PDB ID {pdb_id} is shorter than expected")
        return None
    
    # Construct the path
    pdb_path = os.path.join(pdb_dir, middle_two, pdb_id_lower, f"{pdb_id_lower}_final.pdb")
    
    if os.path.isfile(pdb_path):
        return pdb_path
    else:
        return None

for i, pdb_id in enumerate(pdb_ids):
    if i % 100 == 0:
        logging.info(f"Processing PDB {i}/{len(pdb_ids)}: {pdb_id}")
    
    # Use the new path resolution function
    pdb_path = get_pdb_path(pdb_id)
    
    if pdb_path is None:
        logging.warning(f"PDB file not found for {pdb_id}")
        continue

    found_pdbs += 1

    try:
        structure = parser.get_structure(pdb_id, pdb_path)
        model = structure[0]

        # Get all water atoms with detailed information
        water_atoms = []
        for atom in model.get_atoms():
            if atom.get_parent().get_resname() in ("HOH", "WAT"):
                res = atom.get_parent()
                chain_id = res.get_parent().id
                water_info = {
                    "pdb_id": pdb_id,
                    "chain": chain_id,
                    "resname": res.get_resname(),
                    "resid": res.id[1],
                    "x": float(atom.coord[0]),
                    "y": float(atom.coord[1]),
                    "z": float(atom.coord[2]),
                    "atom_name": atom.get_name(),
                    "bfactor": float(atom.get_bfactor()),
                    "occupancy": float(atom.get_occupancy())
                }
                all_waters_detailed.append(water_info)
                water_atoms.append(atom)

        # Get binding residues for this PDB
        binding_res = binding_df[binding_df["pdb_id"] == pdb_id]
        if binding_res.empty or len(water_atoms) == 0:
            continue

        # Group binding residues by chain to analyze per pocket
        chains_with_binding = binding_res["chain"].unique()
        
        for chain_id in chains_with_binding:
            chain_binding_res = binding_res[binding_res["chain"] == chain_id]
            
            # Get binding atoms for this chain
            chain_res_map = {(row["chain"], row["resi"]) for _, row in chain_binding_res.iterrows()}
            binding_atoms = [
                atom for chain in model
                for res in chain
                for atom in res
                if (chain.id, str(res.id[1])) in chain_res_map and chain.id == chain_id
            ]

            if not binding_atoms:
                continue

            # Calculate water metrics for this pocket (chain) - ALL WATERS
            ns = NeighborSearch(water_atoms)
            
            # 1) All waters within 3.2 Å (allowing double counting)
            all_close_waters = []
            for atom in binding_atoms:
                close = ns.search(atom.coord, 3.2, level="A")
                all_close_waters.extend(close)
            
            total_waters_double_counted = len(all_close_waters)
            
            # 3) Unique waters within 3.2 Å (no double counting)
            unique_water_ids = set((a.get_parent().id, a.get_parent().get_parent().id) for a in all_close_waters)
            unique_waters_count = len(unique_water_ids)
            
            # 2) Normalized waters per binding site residue
            num_binding_residues = len(chain_binding_res)
            normalized_waters_per_residue = unique_waters_count / num_binding_residues if num_binding_residues > 0 else 0
            
            pocket_water_metrics.append({
                "pdb_id": pdb_id,
                "chain": chain_id,
                "pocket_id": f"{pdb_id}_{chain_id}",
                "num_binding_residues": num_binding_residues,
                "total_waters_double_counted": total_waters_double_counted,
                "normalized_waters_per_residue": normalized_waters_per_residue,
                "unique_waters_count": unique_waters_count,
                "data_type": "all_waters"
            })
            
            # Now calculate same metrics for EDIAm-filtered waters
            # Note: EDIAm data might be uppercase, so check both cases
            ediam_pdb_waters = ediam_df[
                (ediam_df["pdb_id"] == pdb_id.upper()) | 
                (ediam_df["pdb_id"] == pdb_id.lower())
            ]
            
            if not ediam_pdb_waters.empty:
                # Match EDIAm waters to PDB structure waters by chain, resid, and insertion code
                ediam_matched_atoms = []
                
                for _, ediam_row in ediam_pdb_waters.iterrows():
                    # Match by chain, residue number, and insertion code
                    for water_atom in water_atoms:
                        water_res = water_atom.get_parent()
                        water_chain = water_res.get_parent().id
                        water_resid = water_res.id[1]
                        water_icode = water_res.id[2].strip() if water_res.id[2].strip() else ""
                        
                        if (water_chain == ediam_row['chain'] and 
                            water_resid == ediam_row['resid'] and
                            water_icode == ediam_row['insertion_code']):
                            ediam_matched_atoms.append(water_atom)
                            break

                if ediam_matched_atoms:
                    # Calculate distances for EDIAm waters
                    ediam_close_waters = []
                    ediam_unique_waters = set()
                    
                    for binding_atom in binding_atoms:
                        for ediam_atom in ediam_matched_atoms:
                            dist = np.linalg.norm(binding_atom.coord - ediam_atom.coord)
                            if dist <= 3.2:
                                ediam_close_waters.append(ediam_atom)
                                unique_id = (ediam_atom.get_parent().id, ediam_atom.get_parent().get_parent().id)
                                ediam_unique_waters.add(unique_id)
                    
                    ediam_total_double = len(ediam_close_waters)
                    ediam_unique_count = len(ediam_unique_waters)
                    ediam_normalized = ediam_unique_count / num_binding_residues if num_binding_residues > 0 else 0
                    
                    pocket_water_metrics.append({
                        "pdb_id": pdb_id,
                        "chain": chain_id,
                        "pocket_id": f"{pdb_id}_{chain_id}",
                        "num_binding_residues": num_binding_residues,
                        "total_waters_double_counted": ediam_total_double,
                        "normalized_waters_per_residue": ediam_normalized,
                        "unique_waters_count": ediam_unique_count,
                        "data_type": "ediam_filtered"
                    })

    except Exception as e:
        logging.error(f"Error processing {pdb_id}: {e}")

logging.info(f"Found valid PDBs for {found_pdbs} / {len(pdb_ids)} structures.")

# === SAVE CSVs ===
logging.info("Saving detailed water data...")
detailed_waters_df = pd.DataFrame(all_waters_detailed)
detailed_waters_df.to_csv(os.path.join(output_dir, "all_waters_detailed.csv"), index=False)

logging.info("Saving pocket water metrics...")
pocket_metrics_df = pd.DataFrame(pocket_water_metrics)
pocket_metrics_df.to_csv(os.path.join(output_dir, "pocket_water_metrics.csv"), index=False)

# === MERGE WITH RIGIDITY DATA ===
logging.info("Merging with rigidity data...")
# Merge on both pdb_id and chain since rigidity data is per chain per pdb
merged_pocket = pd.merge(pocket_metrics_df, rigid_df, on=["pdb_id", "chain"], how="inner")

# === RIGIDITY ORDER ===
cat_order = ["very flexible", "flexible", "rigid", "very rigid"]
rigidity_dtype = CategoricalDtype(categories=cat_order, ordered=True)
merged_pocket["rigidity_quartile"] = merged_pocket["rigidity_quartile"].astype(rigidity_dtype)

# === USE ACTUAL REVERSED VIRIDIS COLORMAP ===
viridis_r_palette = sns.color_palette("viridis_r", 4)

def plot_box_with_kruskal(df, ycol, filename, ylabel, data_type_filter=None):
    """Plot boxplot and perform Kruskal-Wallis test"""
    data = df.copy()
    if data_type_filter:
        data = data[data["data_type"] == data_type_filter]
    
    if data.empty:
        logging.warning(f"No data for {filename} with filter {data_type_filter}")
        return
    
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=data, x="rigidity_quartile", y=ycol, palette=viridis_r_palette, order=cat_order)
    plt.xlabel("")  # no x-axis label
    plt.xticks(ticks=range(4), labels=[label.title() for label in cat_order], fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel(ylabel, fontsize=18)
    
    # Add sample size info for EDIAm plots
    if data_type_filter == "ediam_filtered":
        plt.title(f"EDIAm > 0.3 Waters (n={len(data)} pockets from {data['pdb_id'].nunique()} PDBs)", fontsize=14)
    
    plt.tight_layout()
    
    suffix = f"_{data_type_filter}" if data_type_filter else ""
    plt.savefig(os.path.join(output_dir, f"{filename}{suffix}.png"), dpi=300)
    plt.close()

    # Kruskal-Wallis test
    groups = [grp[ycol].dropna().values for _, grp in data.groupby("rigidity_quartile")]
    groups = [g for g in groups if len(g) > 0]  # Remove empty groups
    
    if len(groups) >= 2:
        stat, p = kruskal(*groups)
        logging.info(f"Kruskal–Wallis for {filename}{suffix}: H={stat:.3f}, p={p:.3e} (n={len(data)} pockets)")
        
        # Save test results
        test_results = pd.DataFrame([{
            "metric": filename,
            "data_type": data_type_filter or "all",
            "kruskal_h": stat,
            "p_value": p,
            "significant": p < 0.05,
            "n_pockets": len(data),
            "n_pdbs": data['pdb_id'].nunique() if 'pdb_id' in data.columns else "N/A"
        }])
        test_file = os.path.join(output_dir, f"kruskal_results_{filename}{suffix}.csv")
        test_results.to_csv(test_file, index=False)
    else:
        logging.info(f"Insufficient groups for Kruskal–Wallis: {filename}{suffix}")

# === PLOT ALL METRICS FOR BOTH ALL WATERS AND EDIAM FILTERED ===
logging.info("Creating plots and statistical tests...")
metrics = [
    ("total_waters_double_counted", "boxplot_total_waters_double_counted", "Total Waters (Double Counted)"),
    ("normalized_waters_per_residue", "boxplot_normalized_waters_per_residue", "Waters per Binding Site Residue"), 
    ("unique_waters_count", "boxplot_unique_waters_count", "Unique Waters Near Binding Site")
]

for col, filename, ylabel in metrics:
    # All waters (main figures)
    plot_box_with_kruskal(merged_pocket, col, filename, ylabel, "all_waters")
    # EDIAm filtered waters (supplementary figures)
    plot_box_with_kruskal(merged_pocket, col, f"supplementary_{filename}", ylabel, "ediam_filtered")

# === SUMMARY STATISTICS ===
logging.info("=== SUMMARY STATISTICS ===")
logging.info(f"Total waters detailed: {len(detailed_waters_df)}")
logging.info(f"Unique PDBs with waters: {detailed_waters_df['pdb_id'].nunique()}")

for data_type in ["all_waters", "ediam_filtered"]:
    subset = merged_pocket[merged_pocket["data_type"] == data_type]
    if not subset.empty:
        logging.info(f"\n{data_type.upper()} ANALYSIS:")
        logging.info(f"Total pockets analyzed: {len(subset)}")
        for col, _, ylabel in metrics:
            mean_val = subset[col].mean()
            std_val = subset[col].std()
            logging.info(f"{ylabel}: mean={mean_val:.2f}, std={std_val:.2f}")

logging.info("Complete water analysis script finished.")
