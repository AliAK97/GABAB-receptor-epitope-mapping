#!/bin/bash
set -euo pipefail

# ============================================================
# Pipeline Configuration — fill in before running
# ============================================================
STRUCTURE_NAME="gaba124_no_temp"

# Step 2 — file name in sorted_output/experimental_structure/
EXP_STRUCTURE="6UO8_chainB.cif"

# Steps 3-5 configuration
ANTIGEN_CHAIN="B"
UNIPROT_CORRECTION=41

ANTIGEN_CHAIN_2=""
UNIPROT_CORRECTION_2=0

THRESHOLD_ALL_PREDS=0.80 
THRESHOLD_TOP_PREDS=0.80
TOP_N=30

# Step 5 (alanine scan)
CHUNK_SIZE=3
TRIALS=5
RELAX_ITER=5
RELAX_NSTRUCT=10
CUTOFF=12
SEED=42

# Step 6 (filter and normalize)
PLDDT_CUTOFF=50
CONTACT_FRACTION=0.10
SALT_BRIDGE_FRACTION=0.10

# ============================================================
# Step Toggles — 1 to run, 0 to skip
# ============================================================
RUN_SEPARATION_RENAME=0      # Step 1 — separate and rename AF3 output files
RUN_DATA_EXTRACTOR=0         # Step 2 — run data extraction
RUN_EPITOPE_FINDER=0         # Step 3 — find common epitope residues
RUN_SALT_BONDS=0             # Step 4 — find salt bond frequencies
RUN_ALANINE_SCAN=0           # Step 5 — alanine mutation scan
RUN_FILTER_NORMALIZE=0       # Step 6 — pLDDT, filter and normalize the epitope table

# ============================================================
# Paths
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${HOME}/.conda/envs/ab_pred_venv/bin/python"
CONDA_ENV="ab_pred_venv"

if [ -f /etc/profile.d/conda.sh ]; then
    set +u
    source /etc/profile.d/conda.sh
    conda activate "$CONDA_ENV"
    set -u
    echo "conda env: ${CONDA_PREFIX:-<none>}"
else
    echo "WARNING: /etc/profile.d/conda.sh not found — running without activation."
    echo "         Steps needing hmmscan (CDR labelling) may be skipped."
fi

# ============================================================
# Validation
# ============================================================
if [ -z "$STRUCTURE_NAME" ]; then
    echo "Error: STRUCTURE_NAME is not set."
    exit 1
fi

REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../" && pwd)"
if [[ "$STRUCTURE_NAME" != /* ]]; then
    STRUCTURE_NAME="${REPO_ROOT}/${STRUCTURE_NAME}"
fi

NAME="$(basename "$STRUCTURE_NAME")"

mkdir -p "${STRUCTURE_NAME}/full_antigen/error_out"

AG2_ARGS=()
if [ -n "$ANTIGEN_CHAIN_2" ]; then
    AG2_ARGS=(--antigen_chain_2 "$ANTIGEN_CHAIN_2"
              --uniprot_correction_2 "$UNIPROT_CORRECTION_2")
fi

echo "Structure: ${NAME}"
if [ -n "$ANTIGEN_CHAIN_2" ]; then
    echo "Two antigen chains: ${ANTIGEN_CHAIN} + ${ANTIGEN_CHAIN_2}"
else
    echo "Single antigen chain: ${ANTIGEN_CHAIN}"
fi

# ============================================================
# Step 1 — separate and rename AF3 output files
# ============================================================
if [[ $RUN_SEPARATION_RENAME -eq 1 ]]; then
    echo "=== Step 1: 01_separation_rename_pdbmaker.py ==="
    "$PYTHON" "${SCRIPT_DIR}/01_separation_rename_pdbmaker.py" "$STRUCTURE_NAME"
fi

# ============================================================
# Step 2 — run data extraction
# ============================================================
if [[ $RUN_DATA_EXTRACTOR -eq 1 ]]; then
    echo "=== Step 2: 02_data_extractor.sh (SLURM job) ==="
    sbatch \
        --output="${STRUCTURE_NAME}/full_antigen/error_out/data_extractor_%j.log" \
        --error="${STRUCTURE_NAME}/full_antigen/error_out/data_extractor_%j.log" \
        "${SCRIPT_DIR}/02_data_extractor.sh" "$STRUCTURE_NAME" "$SCRIPT_DIR" "$EXP_STRUCTURE"
fi

# ============================================================
# Step 3 — find common epitope residues
# ============================================================
if [[ $RUN_EPITOPE_FINDER -eq 1 ]]; then
    echo "=== Step 3: 03_find_common_epitope_residues.py ==="
    "$PYTHON" "${SCRIPT_DIR}/03_find_common_epitope_residues.py" \
        "$STRUCTURE_NAME" \
        --antigen_chain         "$ANTIGEN_CHAIN" \
        --uniprot_correction    "$UNIPROT_CORRECTION" \
        ${AG2_ARGS[@]+"${AG2_ARGS[@]}"} \
        --threshold_all_preds   "$THRESHOLD_ALL_PREDS" \
        --threshold_top_preds   "$THRESHOLD_TOP_PREDS" \
        --top_n                 "$TOP_N"
fi

# ============================================================
# Step 4 — find salt bond frequencies
# ============================================================
if [[ $RUN_SALT_BONDS -eq 1 ]]; then
    echo "=== Step 4: 04_find_salt_bonds.py ==="
    "$PYTHON" "${SCRIPT_DIR}/04_find_salt_bonds.py" \
        "$STRUCTURE_NAME" \
        --antigen_chain        "$ANTIGEN_CHAIN" \
        --uniprot_correction   "$UNIPROT_CORRECTION" \
        ${AG2_ARGS[@]+"${AG2_ARGS[@]}"} \
        --top_n                "$TOP_N" \
        --partner_fraction     "$SALT_BRIDGE_FRACTION"
fi

# ============================================================
# Step 5 — alanine mutation scan
# ============================================================
if [[ $RUN_ALANINE_SCAN -eq 1 ]]; then
    echo "=== Step 5: 05_alanine_scan.sh (SLURM array job) ==="
    LOG="${STRUCTURE_NAME}/full_antigen/error_out/alanine_scan_${NAME}_%j.log"

    ARRAY_MAX=$(( (TOP_N + CHUNK_SIZE - 1) / CHUNK_SIZE - 1 ))
    echo "    ${TOP_N} structures / ${CHUNK_SIZE} per task -> --array=0-${ARRAY_MAX}"

    sbatch \
        --output="$LOG" \
        --error="$LOG" \
        --array="0-${ARRAY_MAX}" \
        --export="ALL,TOP_N=${TOP_N},CHUNK_SIZE=${CHUNK_SIZE},TRIALS=${TRIALS},RELAX_ITER=${RELAX_ITER},RELAX_NSTRUCT=${RELAX_NSTRUCT},CUTOFF=${CUTOFF},SEED=${SEED}" \
        "${SCRIPT_DIR}/05_alanine_scan.sh" \
        "$STRUCTURE_NAME" "$SCRIPT_DIR" "$ANTIGEN_CHAIN" "$UNIPROT_CORRECTION" \
        "$ANTIGEN_CHAIN_2" "$UNIPROT_CORRECTION_2"
fi

# ============================================================
# Step 6 — pLDDT, filter and normalize the epitope table
# ============================================================
if [[ $RUN_FILTER_NORMALIZE -eq 1 ]]; then
    echo "=== Step 6: 06_filter_and_normalize.py ==="
    "$PYTHON" "${SCRIPT_DIR}/06_filter_and_normalize.py" \
        "$STRUCTURE_NAME" \
        --antigen_chain        "$ANTIGEN_CHAIN" \
        --uniprot_correction   "$UNIPROT_CORRECTION" \
        ${AG2_ARGS[@]+"${AG2_ARGS[@]}"} \
        --top_n                "$TOP_N" \
        --plddt_cutoff         "$PLDDT_CUTOFF" \
        --contact_fraction     "$CONTACT_FRACTION" \
        --salt_bridge_fraction "$SALT_BRIDGE_FRACTION"
fi

echo "=== Done ==="
