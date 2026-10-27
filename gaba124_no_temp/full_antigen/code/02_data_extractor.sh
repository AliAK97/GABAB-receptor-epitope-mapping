#!/bin/bash
#SBATCH --job-name=data_extraction_full
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Usage: $0 <structure_name> <code_dir> <exp_structure>"
    exit 1
fi

STRUCTURE_NAME="$1"
CODE_DIR="$2"
EXP_STRUCTURE="$3"

source /etc/profile.d/conda.sh
conda activate ab_pred_venv

~/.conda/envs/ab_pred_venv/bin/python "${CODE_DIR}/02_data_extraction_functions.py" "$STRUCTURE_NAME" \
    --exp_structure "$EXP_STRUCTURE"

conda deactivate
