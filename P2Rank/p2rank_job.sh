#!/bin/bash
#SBATCH --job-name=p2rank_job           # Name of the job
#SBATCH --output=p2rank_job.out        # Output file
#SBATCH --error=p2rank_job.err          # Error file
#SBATCH --time=16:00:00                 # Time limit (adjust as needed)
#SBATCH --mem=25G                       # Memory allocation (adjust as needed)
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --array=1-1%100                 # Set up a job array

# Source necessary environments
source /dors/wankowicz_lab/phenix-installer-dev-5366-intel-linux-2.6-x86_64-centos6/phenix-dev-5366/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

source /dors/wankowicz_lab/shared/conda/etc/profile.d/conda.sh
conda activate qfit

# Define the base directory where your PDB files are stored on the cluster
PDB_DIR="/dors/wankowicz_lab/all_pdb/1_10000"

# Define the directory where results will be saved
RESULTS_DIR="/dors/wankowicz_lab/ellas/p2rank_results/1_10000_results"

# Loop through each PDB ID in the pdb_1_10000.txt file
while read pdb_id; do
    # Check if the pdb_id is not empty
    if [ -n "$pdb_id" ]; then
        # Change directory to the folder for the current pdb_id
        cd "${PDB_DIR}/${pdb_id}" || continue

        # Check if the PDB file exists in the directory
        if [ -f "${pdb_id}.pdb" ]; then
            echo "Processing ${pdb_id}.pdb"

            # Run Phenix to remove heteroatoms
            phenix.pdbtools ${pdb_id}.pdb remove='hetero' output.file_name=${pdb_id}_clean.pdb

            # Run P2Rank for the cleaned PDB file and save output to the results directory
            /dors/wankowicz_lab/ellas/p2rank_2.5/prank predict -f "${pdb_id}_clean.pdb" -o "${RESULTS_DIR}"
        else
            echo "PDB file ${pdb_id}.pdb not found in ${PDB_DIR}/${pdb_id}, skipping..."
        fi
    fi
done < /dors/wankowicz_lab/all_pdb/1_10000/pdb_1_10000.txt
