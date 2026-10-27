"""
Builds the final epitope table. Keeps antigen residues that are in contact with
the antibody at 4Å in more than 10% of all predictions and have a mean pLDDT
above 50 across all predictions. Finally, turns contact and salt bridge counts
into fractions, reporting salt bridges only for residues that form one at 4Å in
more than 10% of all predictions. Adds the mean alanine scan ΔΔG (REU) of each
residue as it is, without normalization.

Reads from analysis_results/:
  <name>_residue_contact_counts.csv      (Step 3)
  salt_bridge_frequencies_<name>.csv     (Step 4)
  alanine_scan_<name>_<chain>_chunk*.csv (Step 5)
Writes into analysis_results/:
  <name>_residue_plddt.csv
  <name>_filtered_and_normalized.csv
"""


import os
import glob
import json
import argparse
import pandas as pd
from Bio.PDB import MMCIFParser


RANK_COLUMN = "iptm"
CUTOFF_TAGS = ["4a", "3a", "2.5a"]
KEY_COLUMNS = ["chain", "residue_id", "residue_num", "residue_num_uniprot"]

AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def residue_plddts(cif_path, json_path, chain_ids):
    """Mean atom pLDDT of every residue on the given chain of one prediction."""
    with open(json_path) as f:
        data = json.load(f)
    atom_plddts = data["atom_plddts"]
    atom_chain_ids = data["atom_chain_ids"]

    model = MMCIFParser(QUIET=True, auth_residues=False).get_structure("mol", cif_path)[0]
    atoms = list(model.get_atoms())

    if len(atoms) != len(atom_plddts) or any(
            atom.get_parent().get_parent().id != chain for atom, chain in zip(atoms, atom_chain_ids)):
        raise ValueError(f"Atoms in {cif_path} do not line up with {json_path}.")

    values = {}
    for atom, plddt in zip(atoms, atom_plddts):
        res = atom.get_parent()
        chain_id = res.get_parent().id
        if chain_id in chain_ids:
            values.setdefault((chain_id, res.id[1], res.resname.strip()), []).append(plddt)
    return {key: sum(v) / len(v) for key, v in values.items()}


def plddt_table(df, cif_dir, pae_dir, antigen_chains):
    """Mean pLDDT per antigen residue across all predictions."""
    chain_ids = [chain for chain, _ in antigen_chains]
    per_prediction = {
        sample: residue_plddts(os.path.join(cif_dir, f"{sample}.cif"),
                               os.path.join(pae_dir, f"{sample}_pae.json"), chain_ids)
        for sample in df["sample_name"]
    }
    plddt_df = pd.DataFrame(per_prediction)
    
    table = plddt_df.mean(axis=1).rename("plddt_mean_all").to_frame()
    table.index.names = ["chain", "residue_num", "residue_id"]
    table = table.reset_index().sort_values(["chain", "residue_num"])
    table["residue_num_uniprot"] = table["residue_num"] + table["chain"].map(dict(antigen_chains))
    return table[KEY_COLUMNS + ["plddt_mean_all"]].round(2)


def ddg_table(analysis_dir, name, antigen_chains, ddg_column):
    """Mean alanine scan ΔΔG (REU) per antigen residue, from every Step 5 chunk."""
    rows = []
    for chain, _ in antigen_chains:
        paths = sorted(glob.glob(os.path.join(analysis_dir, f"alanine_scan_{name}_{chain}_chunk*.csv"))
                       + glob.glob(os.path.join(analysis_dir, f"alanine_scan_{name}_{chain}.csv")))
        if not paths:
            raise SystemExit(f"ERROR: no alanine scan results for chain {chain} in {analysis_dir}")
        scans = pd.concat([pd.read_csv(path, index_col="source_pdb") for path in paths])
        print(f"Chain {chain}: ΔΔG from {len(paths)} files, {len(scans)} structures")

        # A residue is only scanned in structures where it lies within the interface
        # cutoff, so its mean covers those structures only; mean() skips the blanks.
        for mutation, values in scans.items():
            rows.append({
                "chain": chain,
                "residue_num_uniprot": int(mutation[1:-1]),   ### 'K66A' -> 66
                ddg_column: round(values.mean(), 2),
            })
    return pd.DataFrame(rows)


def filter_and_normalize(counts, salt_bridges, plddt, ddg, n_all, top_n, plddt_cutoff,
                         contact_fraction, salt_bridge_fraction):
    """
    Keeps residues with a mean pLDDT across all predictions above 'plddt_cutoff'
    and in contact at 4Å cutoff in more than 'contact_fraction' of all predictions.
    Contact and salt bridge counts become fractions of all predictions or the top N.
    Salt bridge cells stay blank unless the residue forms one at 4Å in more than
    'salt_bridge_fraction' of all predictions. ΔΔG is blank for unscanned residues.
    """
    table = (counts.merge(plddt, on=KEY_COLUMNS, how="left")
                   .merge(salt_bridges, on=KEY_COLUMNS, how="left")
                   .merge(ddg, on=["chain", "residue_num_uniprot"], how="left"))
    keep = ((table["plddt_mean_all"] > plddt_cutoff) &
            (table["count_4a_all"] > contact_fraction * n_all))
    table = table[keep].copy()

    top_suffix = f"top{top_n}_{RANK_COLUMN}"
    fraction_columns = []
    for prefix, new_prefix in (("count", "contact"), ("salt_bridges", "salt_bridge")):
        for suffix, total in (("all", n_all), (top_suffix, top_n)):
            for tag in CUTOFF_TAGS:
                column = f"{new_prefix}_{tag}_{suffix}"
                table[column] = (table[f"{prefix}_{tag}_{suffix}"].fillna(0) / total).round(4)
                fraction_columns.append(column)

    table["residue"] = (table["residue_id"].map(AA3_TO_1).fillna(table["residue_id"])
                        + table["residue_num_uniprot"].astype(str))

    # Partners without their counts, e.g. 'D59(H):503; D1(L):12' -> 'D59(H), D1(L)'.
    table["salt_bridge_partners_4a_all"] = table["partners_4a_all"].fillna("").map(
        lambda partners: ", ".join(p.split(":")[0] for p in partners.split("; ") if p))

    # Blank, not 0, so rare salt bridges stay uncolored in the heatmap.
    rare = table["salt_bridge_4a_all"] <= salt_bridge_fraction
    salt_columns = [c for c in fraction_columns if c.startswith("salt_bridge_")]
    table.loc[rare, salt_columns] = float("nan")
    table.loc[rare, "salt_bridge_partners_4a_all"] = ""

    # Same order as the heatmap: contacts, ΔΔG, then salt bridges.
    ddg_column = [c for c in ddg.columns if c.startswith("ddg_")]
    contact_columns = [c for c in fraction_columns if c.startswith("contact_")]
    return table[["chain", "residue", "plddt_mean_all"] + contact_columns + ddg_column
                 + salt_columns + ["salt_bridge_partners_4a_all"]]


if __name__ == "__main__":
    # ── 1. Arguments and antigen chains ─────────────
    parser = argparse.ArgumentParser(description="Filter and normalize the epitope residue table")
    parser.add_argument("structure_name",
                        help="Path to the structure directory")
    parser.add_argument("--antigen_chain", type=str, required=True,
                        help="Chain ID of the antigen in the CIF files")
    parser.add_argument("--uniprot_correction", type=int, default=0,
                        help="Integer offset added to residue IDs for UniProt numbering")
    parser.add_argument("--antigen_chain_2", type=str, default="",
                        help="Second antigen chain for dimer runs. Empty for a monomer.")
    parser.add_argument("--uniprot_correction_2", type=int, default=0,
                        help="UniProt offset for --antigen_chain_2")
    parser.add_argument("--top_n", type=int, default=30,
                        help="Number of top predictions by iptm, as in Steps 3 and 4")
    parser.add_argument("--plddt_cutoff", type=float, default=50,
                        help="Min mean pLDDT over all predictions (0-100 scale)")
    parser.add_argument("--contact_fraction", type=float, default=0.10,
                        help="Min fraction of all predictions with a 4 Å contact")
    parser.add_argument("--salt_bridge_fraction", type=float, default=0.10,
                        help="Min fraction of all predictions with a 4 Å salt bridge to report it")
    args = parser.parse_args()

    antigen_chains = [(args.antigen_chain, args.uniprot_correction)]
    if args.antigen_chain_2:
        antigen_chains.append((args.antigen_chain_2, args.uniprot_correction_2))
    print(f"Antigen chains: {antigen_chains}")

    # ── 2. Paths ─────────────
    structure_name = str(args.structure_name).rstrip("/")
    name = os.path.basename(structure_name)
    sorted_output = f'{structure_name}/full_antigen/sorted_output'
    csv_path      = f'{sorted_output}/csv_files_metrics/{name}.csv'
    cif_dir       = f'{sorted_output}/cif_folder/{name}_cif'
    pae_dir       = f'{sorted_output}/pae_folder/{name}_json_pae'
    analysis_dir  = f'{structure_name}/full_antigen/analysis_results'
    counts_path   = os.path.join(analysis_dir, f"{name}_residue_contact_counts.csv")
    salt_path     = os.path.join(analysis_dir, f"salt_bridge_frequencies_{name}.csv")

    for path in (csv_path, counts_path, salt_path):
        if not os.path.exists(path):
            raise SystemExit(f"ERROR: input not found: {path}")

    print(f"\n=== Processing structure: {name} ===")
    df = pd.read_csv(csv_path)

    # ── 3. Mean pLDDT per residue ─────────────
    plddt = plddt_table(df, cif_dir, pae_dir, antigen_chains)
    plddt_path = os.path.join(analysis_dir, f"{name}_residue_plddt.csv")
    plddt.to_csv(plddt_path, index=False)
    print(f"Residue pLDDT written to: {plddt_path}")

    # ── 4. Mean ΔΔG per residue ─────────────
    ddg = ddg_table(analysis_dir, name, antigen_chains,
                    ddg_column=f"ddg_mean_top{args.top_n}_{RANK_COLUMN}")

    # ── 5. Filter and normalize ─────────────
    table = filter_and_normalize(
        pd.read_csv(counts_path), pd.read_csv(salt_path), plddt, ddg,
        n_all=len(df), top_n=args.top_n,
        plddt_cutoff=args.plddt_cutoff, contact_fraction=args.contact_fraction,
        salt_bridge_fraction=args.salt_bridge_fraction)

    table_path = os.path.join(analysis_dir, f"{name}_filtered_and_normalized.csv")
    table.to_csv(table_path, index=False)
    print(f"{len(table)} residues kept: {table['residue'].tolist()}")
    print(f"Filtered and normalized table written to: {table_path}")
