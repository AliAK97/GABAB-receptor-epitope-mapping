"""
Find the frequency of epitope residues across the predictions of one structure.

Reads the CSV file written by 02_data_extraction_functions.py and writes,
per antigen chain the frequency of antigen residues closer to antibody 
residues in all and top iptm ranked predictions.
"""

import ast
import argparse
import os
import sys
from collections import Counter
import pandas as pd


CONTACT_CUTOFF_TAGS = ["4a", "3a", "2.5a"]


class _Tee:
    """Writes everything printed to both the terminal and the log file."""
    def __init__(self, terminal, log_file):
        self.terminal = terminal
        self.log_file = log_file

    def write(self, msg):
        """Writes msg to both outputs."""
        self.terminal.write(msg)
        self.log_file.write(msg)

    def flush(self):
        """Flushes both outputs."""
        self.terminal.flush()
        self.log_file.flush()


def optional_int(value):
    """
    Parse an int flag that may arrive empty, e.g. --uniprot_correction_2 "" for
    a single-chain antigen. Returns None so the caller can fall back to a default.
    """
    return int(value) if value.strip() else None


def contact_column(suffix, tag):
    """Name of the Step 2 contact column for an antigen chain suffix and cutoff tag."""
    cutoff_tag = "" if tag == "4a" else f"_{tag}"
    return f"contact_residues{suffix}{cutoff_tag}"


def parse_residue_cell(val):
    """Turns text from the CSV file back into a list of dicts for downstream analysis."""
    if isinstance(val, str):
        try:
            val = ast.literal_eval(val)
        except (ValueError, SyntaxError):
            return []
    return val if isinstance(val, list) else []


def find_common_epitope_residues(df, label, column="contact_residues",
                                 threshold_all_preds=0.8, threshold_top_preds=0.8, top_n=30):
    """Prints the common residue counts in all, highest and lowest iptm ranked structures."""
    def count_residues(subset, thr):
        """Residues in contact in more than a fraction 'thr' of the predictions."""
        counts = Counter()
        total = len(subset)
        if total == 0:
            print(f"[{label}] ERROR: no predictions to count in column '{column}' - epitope list will be empty.")
            return []
        for val in subset[column]:
            seen = {(r['res_id'], r['res_name']) for r in parse_residue_cell(val)}
            counts.update(seen)
        return [
            (res_id, res_name)
            for (res_id, res_name), count in counts.items()
            if (count / total) > thr
        ]

    df_highest = df.nlargest(top_n, "iptm")
    # Bottom N by iptm, for comparison only; not written to the CSV.
    df_lowest  = df.nsmallest(top_n, "iptm")

    all_residue_ids = set()

    def top_counts(subset, sublabel, n=10):
        """Prints the n residues most often in contact."""
        counts = Counter()
        total = len(subset)
        for val in subset[column]:
            counts.update({(r['res_id'], r['res_name']) for r in parse_residue_cell(val)})
        print(f"[{label}] Top {n} most frequent residues ({sublabel}):")
        for (res_id, res_name), cnt in counts.most_common(n):
            print(f"  {res_name}:{res_id}  {cnt}/{total} = {cnt/total:.1%}")

    top_counts(df,         "all rows")
    top_counts(df_highest, f"top {top_n} iptm")
    print("----------------------------")

    print(f"[{label}] Residues present in >{threshold_all_preds*100:.0f}% of all rows:")
    for res_id, res_name in count_residues(df, threshold_all_preds):
        print(f"  {res_name}:{res_id}")
        all_residue_ids.add(res_id)

    print("----------------------------")

    print(f"[{label}] Residues in top {top_n} by iptm:")
    for res_id, res_name in count_residues(df_highest, threshold_top_preds):
        print(f"  {res_name}:{res_id}")
        all_residue_ids.add(res_id)

    print(f"[{label}] Residues in bottom {top_n} by iptm:")
    for res_id, res_name in count_residues(df_lowest, threshold_top_preds):
        print(f"  {res_name}:{res_id}")
        all_residue_ids.add(res_id)

    return sorted(all_residue_ids)


def count_residues_by_distance(df, antigen_specs, top_n=30):
    """
    For every antigen residue, the number of predictions with a contact at each
    distance cutoff, over all predictions and over the top N by iptm. Returns one
    row per residue per chain.
    """
    df_top = df.nlargest(top_n, "iptm")

    def tally(subset, suffix):
        """Contact counts per cutoff, and every residue seen in contact."""
        counters = {tag: Counter() for tag in CONTACT_CUTOFF_TAGS}
        seen_residues = set()
        for tag in CONTACT_CUTOFF_TAGS:
            column = contact_column(suffix, tag)
            if column not in subset.columns:
                continue
            for val in subset[column]:
                seen = {(r['res_id'], r['res_name']) for r in parse_residue_cell(val)}
                seen_residues.update(seen)
                counters[tag].update(seen)
        return counters, seen_residues

    rows = []
    for suffix, chain_id, correction in antigen_specs:
        counters_all, residues_all = tally(df, suffix)
        counters_top, residues_top = tally(df_top, suffix)

        for res_id, res_name in sorted(residues_all):
            key = (res_id, res_name)
            row = {
                "chain": chain_id,
                "residue_id": res_name,
                "residue_num": res_id,
                "residue_num_uniprot": res_id + correction,
            }
            for tag in CONTACT_CUTOFF_TAGS:
                row[f"count_{tag}_all"] = counters_all[tag].get(key, 0)
                row[f"count_{tag}_top{top_n}_iptm"] = (
                    counters_top[tag].get(key, 0) if key in residues_top else pd.NA)
            rows.append(row)

    return pd.DataFrame(rows)


if __name__ == "__main__":
    # ── 1. Arguments ─────────────
    parser = argparse.ArgumentParser(
        description="Find common epitope residues for one structure")
    parser.add_argument("structure_name",
                        help="Path to the structure directory")
    parser.add_argument("--antigen_chain", type=str, default="A",
                        help="Chain ID of the first antigen")
    parser.add_argument("--antigen_chain_2", type=str, default=None,
                        help="Chain ID of the second antigen (optional)")
    parser.add_argument("--uniprot_correction", type=int, default=0,
                        help="Integer adjustment to match UniProt numbering (chain 1)")
    parser.add_argument("--uniprot_correction_2", type=optional_int, default=None,
                        help="Integer adjustment to match UniProt numbering (chain 2)")
    parser.add_argument("--threshold_all_preds", type=float, default=0.8,
                        help="Min frequency in all rows to print out a residue")
    parser.add_argument("--threshold_top_preds", type=float, default=0.8,
                        help="Min frequency in top/bottom N print out a residue")
    parser.add_argument("--top_n", type=int, default=30,
                        help="Number of top/bottom predictions by iptm used for both the ranked "
                             "epitope lists and the contact counts")
    args = parser.parse_args()

    correction_2 = (args.uniprot_correction_2 if args.uniprot_correction_2 is not None
                    else args.uniprot_correction)

    # ── 2. Paths and antigen chains ─────────────
    structure_name = str(args.structure_name).rstrip("/")
    name = os.path.basename(structure_name)
    sorted_output = f'{structure_name}/full_antigen/sorted_output'
    csv_path     = f'{sorted_output}/csv_files_metrics/{name}.csv'
    analysis_dir = f'{structure_name}/full_antigen/analysis_results'
    os.makedirs(analysis_dir, exist_ok=True)

    antigen_specs = [("", args.antigen_chain, args.uniprot_correction)]
    if args.antigen_chain_2:
        antigen_specs.append(("_ag2", args.antigen_chain_2, correction_2))

    # Everything printed from here on also goes to epitope_finder_debug.txt.
    log_path = os.path.join(analysis_dir, "epitope_finder_debug.txt")
    log_file = open(log_path, "w")
    sys.stdout = _Tee(sys.__stdout__, log_file)

    try:
        if not os.path.exists(csv_path):
            print(f"Skipping {name}: {csv_path} not found")
            sys.exit(1)

        df = pd.read_csv(csv_path)

        present_specs = []
        for spec in antigen_specs:
            suffix, chain_id, _ = spec
            column = contact_column(suffix, "4a")
            if column in df.columns:
                present_specs.append(spec)
            else:
                print(f"WARNING: chain {chain_id} requested but '{column}' is not in {csv_path}")
        antigen_specs = present_specs

        if not antigen_specs:
            print(f"ERROR: none of the requested antigen chains are present in {csv_path}.")
            sys.exit(1)

        # ── 3. Common epitope residues per antigen chain ─────────────
        for suffix, chain_id, correction in antigen_specs:
            residue_ids = find_common_epitope_residues(
                df,
                label=f"{name}_{chain_id}",
                column=contact_column(suffix, "4a"),
                threshold_all_preds=args.threshold_all_preds,
                threshold_top_preds=args.threshold_top_preds,
                top_n=args.top_n)
            residue_nums = [r + correction for r in residue_ids]

            print("----------------------------")
            print(f"[{name}] Chain {chain_id} residue IDs (raw PDB):           {residue_ids}")
            print(f"[{name}] Chain {chain_id} with uniprot_correction ({correction}): {residue_nums}")

        # ── 4. Contact counts by distance ─────────────
        print(f"[{name}] Counting contact residues by distance (top_n={args.top_n})...")
        counts_df = count_residues_by_distance(df, antigen_specs, top_n=args.top_n)
        counts_path = os.path.join(analysis_dir, f"{name}_residue_contact_counts.csv")
        counts_df.to_csv(counts_path, index=False)
        print(f"Residue contact counts written to: {counts_path}")
    finally:
        sys.stdout = sys.__stdout__
        log_file.close()
        print(f"Debug log written to: {log_path}")
