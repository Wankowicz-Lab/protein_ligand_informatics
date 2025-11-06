#!/bin/bash
#SBATCH --job-name=freesasa_peratom
#SBATCH --output=/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/logs/freesasa_peratom_%A_%a.out
#SBATCH --error=/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/logs/freesasa_peratom_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Exit immediately if any command fails
set -e

# === CONFIG ===
PDB_LIST="/panfs/accrepfs.vampire/data/wankowicz_lab/.../pdb_list.txt"
BASE_DIR="/panfs/accrepfs.vampire/data/wankowicz_lab/PDBRedo/pdb-redo3"
OUT_DIR="/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_per_atom"
LOG_FILE="/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_missing_pdbs.log"
FAILED_NO_PDB="/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_failed_no_pdb_files.txt"
FAILED_YES_PDB="/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/freesasa_failed_yes_pdb_files.txt"

mkdir -p "$OUT_DIR"
mkdir -p "$(dirname "$LOG_FILE")"

# Clear previous failed files lists
> "$FAILED_NO_PDB"
> "$FAILED_YES_PDB"

# === ENVIRONMENT ===
module load StdEnv/2020 gcc/9.3.0 freesasa/2.1.0

# === LOOP OVER PDBS ===
while read -r pdb_lc; do
    pdb_lc=$(echo "$pdb_lc" | tr -d '[:space:]')
    [ -z "$pdb_lc" ] && continue

    mid=${pdb_lc:1:2}
    pdb_path="${BASE_DIR}/${mid}/${pdb_lc}/${pdb_lc}_final.pdb"
    out_file="${OUT_DIR}/${pdb_lc}_sasa_atom.json"

    if [ ! -f "$pdb_path" ]; then
        reason="PDB file not found at ${pdb_path}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Missing: $pdb_lc - $reason" | tee -a "$LOG_FILE"
        echo "$pdb_lc" >> "$FAILED_NO_PDB"
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running FreeSASA on $pdb_lc ..."
    
    # Run FreeSASA - redirect stderr to suppress module messages
    freesasa --format=json --depth=atom "$pdb_path" > "$out_file" 2>/dev/null
    exit_code=$?

    if [ $exit_code -ne 0 ]; then
        reason="FreeSASA failed with exit code $exit_code"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FATAL ERROR: $pdb_lc - $reason" | tee -a "$LOG_FILE"
        echo "$pdb_lc" >> "$FAILED_YES_PDB"
        echo "Job terminated due to FreeSASA failure on $pdb_lc"
        exit 1
    fi

    # Check if output file is empty or too small (less than 100 bytes indicates failure)
    if [ ! -s "$out_file" ] || [ $(stat -c%s "$out_file") -lt 100 ]; then
        reason="FreeSASA produced empty or invalid output file"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FATAL ERROR: $pdb_lc - $reason" | tee -a "$LOG_FILE"
        echo "$pdb_lc" >> "$FAILED_YES_PDB"
        echo "Job terminated due to empty/invalid output for $pdb_lc"
        exit 1
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✔ Done: $pdb_lc"
    
done < "$PDB_LIST"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Summary:"
echo "  - Failed (no PDB): $(wc -l < "$FAILED_NO_PDB") entries in $FAILED_NO_PDB"
echo "  - Failed (PDB exists): $(wc -l < "$FAILED_YES_PDB") entries in $FAILED_YES_PDB"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] All jobs completed successfully!"
