# GABAB-receptor-epitope-mapping

Code for predicting the antibody-antigen complex and epitope mapping of a monoclonal antibody against GABA<sub>B</sub> receptor binding to GABA<sub>B</sub> subunit 2 with AlphaFold3.

The antibody Fab (heavy chain `H`, light chain `L`) is co-folded with the GABA<sub>B2</sub> extracellular domain (GB2-ECD) referred to as chain `B` over 120 seeds × 5 diffusion samples, producing 600 models. Once the predictions are ready, the model extracts multiple features from them. Finally across all models and the top-ranked ones by ipTM, it measures:

- how often each antigen residue contacts the antibody residue (at 4, 3 and 2.5 Å);
- which antigen–antibody residues form potential salt bond;
- the binding ΔΔG of each interface residue mutation, using a PyRosetta alanine scan.

These are merged into a single filtered epitope table. A notebook reproduces the structure figures of the manuscript.

> **Citation:** _TODO: add manuscript reference / DOI_

---

## Repository layout

The repository is laid out exactly as the pipeline expects it at runtime, so it can be run straight after cloning.

```
GABAB-receptor-epitope-mapping/
├── data/
│   └── experimental/                   # 4mqe.cif, 4ms4.cif: crystal structures used by the notebook
├── notebooks/
│   └── visualizing_structure.ipynb     
└── gaba124_no_temp/                    # one folder per structure ("structure directory")
    └── full_antigen/
        ├── input/
        │   └── gaba124_no_temp.json    
        ├── run_af3.sh                  # NOT included: write your own (see step 0)
        ├── output/ 
        ├── code/                       # analysis pipeline (steps 1–6)
        ├── sorted_output/
        │   ├── meta_data.csv           
        │   ├── experimental_structure/ 
        │   ├── cif_folder/
        │   ├── json_folder/
        │   ├── pae_folder/
        │   ├── pdb_folder/
        │   └── csv_files_metrics/
        ├── analysis_results/
        └── error_out/
```

The empty folders hold a `.gitkeep` file to be tracked. The scripts also create the folders even if the `.gitkeep` is missing.

---

## Requirements

- **AlphaFold3**: installed separately; see the [AlphaFold3 repository](https://github.com/google-deepmind/alphafold3).
- **PyRosetta**: used for the alanine scan (step 5); see [PyRosetta licensing](https://www.pyrosetta.org/home/licensing-pyrosetta) and the [download page](https://www.pyrosetta.org/downloads).
- **SLURM**: steps 2 and 5 are submitted with `sbatch`.
- **GPU**: NVIDIA A100 80GB GPU and CUDA Version 12.7.

The scripts expect the environment to be called `ab_pred_venv` and to live at `~/.conda/envs/ab_pred_venv`. They activate it through `/etc/profile.d/conda.sh`. If your setup differs, edit `PYTHON` / `CONDA_ENV` in `run_pipeline.sh` and the `conda activate` lines in `02_data_extractor.sh` and `05_alanine_scan.sh`. You may also need to change the `#SBATCH` partition (`compute`) and time limits for your cluster.

---

## Installation

```bash
git clone https://github.com/AliAK97/GABAB-receptor-epitope-mapping
cd GABAB-receptor-epitope-mapping
conda env create -f environment.yml
conda activate ab_pred_venv
```

PyRosetta is not in `environment.yml`. After obtaining a license, install it into the same environment, e.g.:

```bash
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
```

---

## Usage

All paths below are relative to `gaba124_no_temp/full_antigen/`.

### Step 0: acquire input file and run AlphaFold3

- The input file is not part of this repository and will be made available under a Material Transfer Agreement (MTA) for non-commercial usage through the corresponding author.
- The AlphaFold3 run script is also not part of this repository, because it depends on how AlphaFold3 is installed on your system. Create `run_af3.sh` in `full_antigen/`. It must read `input/gaba124_no_temp.json` and write the predictions in `output/`.

### Steps 1–6: configure and run the pipeline

Open `code/run_pipeline.sh` and set the configuration block:

| Variable | Default | Meaning |
|---|---|---|
| `STRUCTURE_NAME` | `gaba124_no_temp` | Structure directory, relative to the repository root (or an absolute path) |
| `EXP_STRUCTURE` | `6UO8_chainB.cif` | Reference structure in `sorted_output/experimental_structure/` |
| `ANTIGEN_CHAIN` | `B` | Antigen chain ID in the predictions |
| `UNIPROT_CORRECTION` | `41` | Offset to translate the sequence index to UniProt index (UniProt = model + 41) |
| `ANTIGEN_CHAIN_2`, `UNIPROT_CORRECTION_2` | empty, `0` | Optional second antigen chain |
| `THRESHOLD_ALL_PREDS`, `THRESHOLD_TOP_PREDS` | `0.80` | Contact frequency above which step 3 prints a residue as epitope |
| `TOP_N` | `30` | Number of top-ipTM models used for the "top" statistics and the alanine scan |
| `CHUNK_SIZE`, `TRIALS`, `RELAX_ITER`, `RELAX_NSTRUCT`, `CUTOFF`, `SEED` | `3`, `5`, `5`, `10`, `12`, `42` | Alanine scan settings (step 5) |
| `PLDDT_CUTOFF`, `CONTACT_FRACTION`, `SALT_BRIDGE_FRACTION` | `50`, `0.10`, `0.10` | Filters for the final table (step 6) |

Then switch steps on (`1`) or off (`0`) in the `RUN_*` toggles and run:

```bash
bash gaba124_no_temp/full_antigen/code/run_pipeline.sh
```

Steps 2 and 5 are submitted as SLURM jobs. Wait for them to finish before running the steps that depend on them. A typical order is:

0. Step 0, running predictions
1. Step 1
2. Step 2, then wait for the job
3. Steps 3, 4 and 5, then wait for the step 5 array
4. Step 6

---

## Notebook

`notebooks/visualizing_structure.ipynb` reproduces Figures 3A–C and Supplementary Figures 1, 3 and 4. It reads:

- the predicted models from `sorted_output/cif_folder/gaba124_no_temp_cif/` (step 1);
- the metrics table `sorted_output/csv_files_metrics/gaba124_no_temp.csv` (step 2);
- the experimental structures in `data/experimental/`.

Run steps 0–2 before opening the notebook.

> **Data availability:** the 600 predicted models and the metrics table are not stored in this repository because of their size. They can be reproduced by running the input file with AF3

---

## Running another structure

1. Copy `gaba124_no_temp/` to a new folder named after your structure, then empty the step-filled folders.
2. Replace `input/<name>.json` with your AlphaFold3 input.
3. Add a row to `sorted_output/meta_data.csv` with the columns `structure_name, heavy_chain, heavy_chain_seq, light_chain, light_chain_seq, antigen_chain, antigen_chain_seq, antigen_chain_2, antigen_chain_2_seq`.
4. Put the reference structure in `sorted_output/experimental_structure/`.
5. Update the variables in `code/run_pipeline.sh`.

Steps 4 and 5 assume the antibody chains are named `H` and `L`.

---

## License

MIT; see [LICENSE](LICENSE).
