#!/bin/bash
#SBATCH --job-name=a_scan_full
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
# NOTE: --array is set by run_pipeline.sh, which derives it from
#       TOP_N and CHUNK_SIZE. Do not hardcode it here.

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    echo "Usage: $0 <structure_name> <code_dir> <antigen_chain> <uniprot_correction>" \
         "[antigen_chain_2] [uniprot_correction_2]"
    exit 1
fi

STRUCTURE_NAME="$1"
CODE_DIR="$2"
ANTIGEN_CHAIN="$3"
UNIPROT_CORRECTION="$4"
# Optional second antigen chain. It is scanned as its own run (partners "B_HL"),
# not merged with the first chain's interface.
ANTIGEN_CHAIN_2="${5:-}"
UNIPROT_CORRECTION_2="${6:-0}"

TOP_N="${TOP_N:-30}"
CHUNK_SIZE="${CHUNK_SIZE:-3}"
CHUNK_INDEX="${SLURM_ARRAY_TASK_ID:-0}"

source /etc/profile.d/conda.sh
conda activate ab_pred_venv

~/.conda/envs/ab_pred_venv/bin/python "${CODE_DIR}/05_alanine_mutation.py" \
    "$STRUCTURE_NAME" \
    --antigen_chain        "$ANTIGEN_CHAIN" \
    --antigen_chain_2      "$ANTIGEN_CHAIN_2" \
    --uniprot_correction   "$UNIPROT_CORRECTION" \
    --uniprot_correction_2 "$UNIPROT_CORRECTION_2" \
    --top_n                "$TOP_N" \
    --trials               "${TRIALS:-5}" \
    --relax_iter           "${RELAX_ITER:-5}" \
    --relax_nstruct        "${RELAX_NSTRUCT:-10}" \
    --cutoff               "${CUTOFF:-12}" \
    --seed                 "${SEED:-42}" \
    --chunk_size           "$CHUNK_SIZE" \
    --chunk_index          "$CHUNK_INDEX"

conda deactivate
