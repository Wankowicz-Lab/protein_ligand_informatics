#!/bin/bash
#SBATCH --mem=48G  # Request 48GB of memory
#SBATCH --cpus-per-task=4  # Request 8 CPUs per task
#SBATCH --time=48:00:00  # Request a walltime of 24 hours
#SBATCH --account=wankowicz_lab
#SBATCH --mail-user=ahmed.j.aslam@vanderbilt.edu
#SBATCH --array=1-25200%200  # Set up a job array
#SBATCH --output=out/proc._%j_stdout
#SBATCH --job-name=lig_clustering_editdistance


module load python/3.8

OUR_DIR="./built_225_largest_adj_mxs"
FILES=($OUR_DIR/*.json)
NUM_FILES=${#FILES[@]}

TASK_ID=$SLURM_ARRAY_TASK_ID

i=$(( (TASK_ID - 1) / (NUM_FILES - 1) ))
j=$(( (TASK_ID - 1) % (NUM_FILES - 1) ))

if [ $j -ge $i ]; then
j=$((j + 1))
fi

python task_get_editdistance.py "${FILES[$i]}" "${FILES[$j]}"