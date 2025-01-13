#!/bin/bash
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00
#SBATCH --account=wankowicz_lab
#SBATCH --mail-user=ahmed.j.aslam@vanderbilt.edu
#SBATCH --array=1-25200%500
#SBATCH --output=out2/proc._%j_stdout
#SBATCH --job-name=lig_clustering_editdistance

module load python/3.8

OUR_DIR="./built_225_largest_adj_mxs"
FILES=($OUR_DIR/*.json)
NUM_FILES=${#FILES[@]}

TASK_ID=$SLURM_ARRAY_TASK_ID

tempIndex=$TASK_ID
i=0
span=$(( (NUM_FILES-1) - i ))

while [ $tempIndex -gt $span ]; do
  tempIndex=$(( tempIndex - span ))
  i=$(( i + 1 ))
  span=$(( (NUM_FILES-1) - i ))
done

j=$(( i + tempIndex ))

echo "TASK_ID=${TASK_ID}: i=${i}, j=${j}"

python task_get_editdistance.py "${FILES[$i]}" "${FILES[$j]}"