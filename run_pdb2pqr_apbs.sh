#!/bin/bash
#SBATCH --job-name=apbs_array
#SBATCH --output=/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/apbs_outputs/logs/slurm_%A_%a.out
#SBATCH --error=/panfs/accrepfs.vampire/data/wankowicz_lab/ellas/apbs_outputs/logs/slurm_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# =========================
# User/configurable paths
# =========================
BASE_DIR="/panfs/accrepfs.vampire/data/wankowicz_lab/ellas"
LIST="${BASE_DIR}/pdb_list.txt"

# Input PDB directories
pdb_dir_primary="/panfs/accrepfs.vampire/data/wankowicz_lab/PDBRedo/pdb-redo3"
pdb_dir_fallback="${BASE_DIR}/pdbredo_pdbs"

# Outputs (ALL intermediates & results go here)
OUT_DIR="${BASE_DIR}/apbs_outputs"
LOG_DIR="${OUT_DIR}/logs"
MAP_FILE="${OUT_DIR}/resname_map.tsv"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# ================
# Environment
# ================
# Ensure APBS is on PATH (installed under ellas/)
export PATH="${BASE_DIR}/apbs/APBS-3.4.1.Linux/bin:${PATH}"

# Activate the conda env that has pdb2pqr
if [ -f "${BASE_DIR}/miniconda3/etc/profile.d/conda.sh" ]; then
  source "${BASE_DIR}/miniconda3/etc/profile.d/conda.sh"
fi
conda activate qfit >/dev/null 2>&1 || true

# Sanity checks
which pdb2pqr >/dev/null 2>&1 || { echo "[FATAL] pdb2pqr not found on PATH"; exit 2; }
which apbs    >/dev/null 2>&1 || { echo "[FATAL] apbs not found on PATH"; exit 2; }

# =========================
# Residue remap table
# =========================
if [ ! -s "${MAP_FILE}" ]; then
  cat > "${MAP_FILE}" <<'EOF'
MSE MET
SEC CYS
HID HIS
HIE HIS
HIP HIS
ASH ASP
GLH GLU
LYN LYS
CYX CYS
CSS CYS
EOF
fi

# =========================
# Helpers
# =========================

# get_pdb_path UPPER
get_pdb_path () {
  local pdb_upper="$1"
  local pdb_lc
  pdb_lc="$(echo "$pdb_upper" | tr '[:upper:]' '[:lower:]')"
  if [ "${#pdb_lc}" -lt 4 ]; then
    echo ""
    return
  fi
  local middle_two="${pdb_lc:1:2}"
  local primary="${pdb_dir_primary}/${middle_two}/${pdb_lc}/${pdb_lc}_final.pdb"
  if [ -s "${primary}" ]; then
    echo "${primary}"
    return
  fi
  local fallback="${pdb_dir_fallback}/${pdb_lc}_final.pdb"
  if [ -s "${fallback}" ]; then
    echo "${fallback}"
    return
  fi
  echo ""
}

# Clean one PDB: remap, drop water/ligands, run pdb2pqr -> apbs
process_one () {
  local pdb_upper="$1"         # e.g., 2DB7
  local pdb_path="$2"          # absolute path to *_final.pdb
  local core_lc
  core_lc="$(echo "$pdb_upper" | tr '[:upper:]' '[:lower:]')"  # '2db7'

  local remap_pdb="${OUT_DIR}/${core_lc}_remap.pdb"
  local clean_pdb="${OUT_DIR}/${core_lc}_clean.pdb"
  local out_pqr="${OUT_DIR}/${core_lc}.pqr"
  local out_in="${OUT_DIR}/${core_lc}.in"
  local log="${LOG_DIR}/${core_lc}.log"
  local err="${LOG_DIR}/${core_lc}.err"
  local timefile="${LOG_DIR}/${core_lc}.time"

  # 1) Remap common variants -> 20 AAs
  awk -v MAP="${MAP_FILE}" '
    BEGIN{
      while((getline<MAP)>0){from=$1; to=$2; map[toupper(from)]=toupper(to)}
      FS=""; OFS=""
    }
    /^(ATOM|HETATM)/{
      res=substr($0,18,3)
      gsub(/^[ ]+|[ ]+$/,"",res)
      up=toupper(res)
      if(up in map){
        $0 = substr($0,1,17) sprintf("%3s", map[up]) substr($0,21)
      }
    }
    {print}
  ' "${pdb_path}" > "${remap_pdb}" 2>>"${err}"

  if [ ! -s "${remap_pdb}" ]; then
    echo "[ERROR] Remapped file empty for ${pdb_upper}" | tee -a "${err}"
    return
  fi

  # 2) Drop waters and ligands; keep simple ions
  awk '
    /^(ATOM)/ {print; next}
    /^(HETATM)/ {
      res=substr($0,18,3); gsub(/^[ ]+|[ ]+$/,"",res)
      if (res=="HOH") next
      if (res ~ /^(NA|CL|ZN|MG|CA)$/) print
      next
    }
  ' "${remap_pdb}" > "${clean_pdb}" 2>>"${err}"

  if [ ! -s "${clean_pdb}" ]; then
    echo "[ERROR] Clean file empty for ${pdb_upper}" | tee -a "${err}"
    return
  fi

  # 3) Run pdb2pqr (PARSE; drop-water; generate .in)
  /usr/bin/time -f "%E %M" -o "${timefile}" \
    pdb2pqr --ff=PARSE --drop-water --apbs-input "${out_in}" \
      "${clean_pdb}" "${out_pqr}" > "${log}" 2>> "${err}"

  if [ ! -s "${out_pqr}" ] || [ ! -s "${out_in}" ]; then
    echo "[ERROR] PDB2PQR failed for ${pdb_upper}" | tee -a "${err}"
    tail -n 50 "${log}" >> "${err}"
    return
  fi

  # 4) Patch .in to use absolute .pqr path
  sed -i "s#^\([[:space:]]*mol[[:space:]]\+pqr[[:space:]]\+\).*#\1${out_pqr}#" "${out_in}"

  # 5) Run APBS
  apbs "${out_in}" >> "${log}" 2>> "${err}"
}

# =========================
# Load list & shard by array index (5 shards)
# =========================
if [ ! -s "${LIST}" ]; then
  echo "[FATAL] PDB list not found: ${LIST}"
  exit 2
fi

mapfile -t PDBS < <(awk 'NF{print toupper($0)}' "${LIST}")

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
SHARDS=5

missing_file="${OUT_DIR}/missing_pdbs_task${TASK_ID}.txt"
: > "${missing_file}"

total="${#PDBS[@]}"
echo "[INFO] Total IDs: ${total} | This task: index mod ${SHARDS} == ${TASK_ID}"

processed=0
found=0
for idx in "${!PDBS[@]}"; do
  # shard selection (covers all IDs even if total % 5 != 0)
  mod=$(( idx % SHARDS ))
  if [ "${mod}" -ne "${TASK_ID}" ]; then
    continue
  fi

  pdb_upper="${PDBS[$idx]}"
  (( processed++ )) || true
  if (( processed % 50 == 0 )); then
    echo "[INFO][task ${TASK_ID}] processed=${processed} (idx=${idx}) last=${pdb_upper}"
  fi

  pdb_path="$(get_pdb_path "${pdb_upper}")"
  if [ -z "${pdb_path}" ]; then
    echo "${pdb_upper}" >> "${missing_file}"
    continue
  fi

  (( found++ )) || true
  process_one "${pdb_upper}" "${pdb_path}"
done

echo "[INFO][task ${TASK_ID}] done. processed_in_shard=${processed}, found_inputs=${found}"
echo "[INFO] Missing list (this shard): ${missing_file}"
