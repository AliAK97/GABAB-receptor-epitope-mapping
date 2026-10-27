"""
This code identifies charged residues on antigen and antibody chains that are closer
to each other than a certain threshold and can potentially form a salt bond.

Writes two files into analysis_results/:
  salt_bridge_frequencies_<name>.txt   every salt bridge found, unfiltered
  salt_bridge_frequencies_<name>.csv   counts from the frequent antibody partners only:
                                       those bridging at 4 Å in more than 10% of all predictions
"""


import pandas as pd
import os
import argparse
from collections import Counter
from Bio.PDB import MMCIFParser
from abnumber import Chain


HC_CHAIN = 'H'
LC_CHAIN = 'L'
CUTOFFS = [4.0, 3.0, 2.5]

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}
CHARGED_ATOMS = {
    'ARG': ['NH1', 'NH2'],
    'LYS': ['NZ'],
    'HIS': ['ND1', 'NE2'],
    'ASP': ['OD1', 'OD2'],
    'GLU': ['OE1', 'OE2'],
}

POSITIVE = {'ARG', 'LYS', 'HIS'}
NEGATIVE = {'ASP', 'GLU'}

RANK_COLUMN = "iptm"

CDR_SCHEME = "imgt"

AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

_PAIR_CACHE = {}


def _charged_residues(chain):
    """Charged residues of a chain, key is the residue id."""
    return {res.id[1]: res for res in chain if res.resname.strip() in CHARGED_ATOMS}


def _min_charged_distance(res_a, res_b):
    """Closest distance between the charged side-chain atoms of two residues in opposite chains."""
    n_a, n_b = res_a.resname.strip(), res_b.resname.strip()
    atoms_a = [res_a[a] for a in CHARGED_ATOMS[n_a] if a in res_a and a not in BACKBONE_ATOMS]
    atoms_b = [res_b[a] for a in CHARGED_ATOMS[n_b] if a in res_b and a not in BACKBONE_ATOMS]
    return min((x - y for x in atoms_a for y in atoms_b), default=None)


def compute_charged_pairs(cif_file, antigen_chain):
    """Returns every oppositely-charged antigen-antibody residue pair with distance."""
    key = (cif_file, antigen_chain)
    if key in _PAIR_CACHE:
        return _PAIR_CACHE[key]

    structure = MMCIFParser(QUIET=True, auth_residues=False).get_structure("mol", cif_file)

    antigen_charged = _charged_residues(structure[0][antigen_chain])
    hc_charged      = _charged_residues(structure[0][HC_CHAIN])
    lc_charged      = _charged_residues(structure[0][LC_CHAIN])

    pairs = []
    for ag_id, ag_res in antigen_charged.items():
        n_ag = ag_res.resname.strip()
        for partner_chain, partner_charged in ((HC_CHAIN, hc_charged), (LC_CHAIN, lc_charged)):
            for p_id, p_res in partner_charged.items():
                n_p = p_res.resname.strip()
                if not ((n_ag in POSITIVE and n_p in NEGATIVE) or
                        (n_ag in NEGATIVE and n_p in POSITIVE)):
                    continue
                dist = _min_charged_distance(ag_res, p_res)
                if dist is not None:
                    pairs.append((partner_chain, ag_id, n_ag, p_id, n_p, dist))

    _PAIR_CACHE[key] = pairs
    return pairs


def find_salt_bridges(cif_file, antigen_chain, cutoff=4.0):
    """
    Returns salt bridges based on the 'cutoff' distance: the antigen, HC and LC
    residues involved, and the antigen-antibody residue pairs themselves.
    """
    antigen_ids, hc_ids, lc_ids, pairs = set(), set(), set(), set()

    for partner_chain, ag_id, n_ag, p_id, n_p, dist in compute_charged_pairs(cif_file, antigen_chain):
        if dist > cutoff:
            continue
        antigen_ids.add((ag_id, n_ag))
        (hc_ids if partner_chain == HC_CHAIN else lc_ids).add((p_id, n_p))
        pairs.add(((ag_id, n_ag), (partner_chain, p_id, n_p)))

    return sorted(antigen_ids), sorted(hc_ids), sorted(lc_ids), sorted(pairs)


def extract_chain_sequence(structure, chain_id):
    """One-letter sequence of a chain, and the residue id at each position."""
    seq, res_ids = [], []
    for res in structure[0][chain_id]:
        aa = AA3_TO_1.get(res.resname.strip())
        if aa is None:
            continue
        seq.append(aa)
        res_ids.append(res.id[1])
    return "".join(seq), res_ids


def derive_cdr_ranges(seq, res_ids, scheme=CDR_SCHEME):
    """CDR ranges as {name: (first_res_id, last_res_id)}, taken from abnumber's numbering."""
    chain = Chain(seq, scheme=scheme)

    domain_offset = seq.find(chain.seq)
    if domain_offset < 0:
        raise ValueError("Variable domain from abnumber not found in the chain sequence.")

    cdr_indices = {'CDR1': [], 'CDR2': [], 'CDR3': []}
    for i, (pos, _) in enumerate(chain):
        region = pos.get_region()
        if region in cdr_indices:
            cdr_indices[region].append(domain_offset + i)

    return {name: (res_ids[idx[0]], res_ids[idx[-1]]) for name, idx in cdr_indices.items()}


def find_cdr_ranges(cif_file, scheme=CDR_SCHEME):
    """CDR ranges of both antibody chains, as {'HC': {name: (first, last)}, 'LC': ...}."""
    structure = MMCIFParser(QUIET=True, auth_residues=False).get_structure("mol", cif_file)

    ranges = {}
    for chain_id, label in ((HC_CHAIN, 'HC'), (LC_CHAIN, 'LC')):
        seq, res_ids = extract_chain_sequence(structure, chain_id)
        ranges[label] = derive_cdr_ranges(seq, res_ids, scheme=scheme)
    return ranges


def format_cdr_ranges(cdr_ranges, scheme=CDR_SCHEME):
    """Report lines listing the CDR ranges, for the log and .txt file."""
    lines = [f"\nCDR ranges ({scheme.upper()} definition, CIF residue numbers):"]
    for chain in ('HC', 'LC'):
        for name, (start, end) in cdr_ranges[chain].items():
            lines.append(f"  {chain} {name}: residues {start}-{end}")
    return lines


def get_cdr_name(res_id, cdr_ranges, chain='HC'):
    """Name of the CDR a residue falls in, or None if it is in FC."""
    prefix = '' if chain == 'HC' else 'LC_'
    for name, (start, end) in cdr_ranges[chain].items():
        if start <= res_id <= end:
            return f"{prefix}{name}"
    return None


def count_salt_bridges(subset, cdr_ranges, cutoff=4.0, label=""):
    """Returns bridge frequencies across predictions for every included antigen chain."""
    n = len(subset)
    bridge_counter = {}
    hc_bridge_counter = {}
    lc_bridge_counter = {}
    residue_counter = Counter()
    partner_counter = {}

    for _, row in subset.iterrows():
        sample = row['sample_name']
        cif    = os.path.join(CIF_DIR, f"{sample}.cif")

        hc_seen, lc_seen = set(), set()

        for antigen_chain, uniprot_correction in ANTIGEN_CHAINS:
            ag_bridges, hc_bridges, lc_bridges, pairs = find_salt_bridges(
                cif, antigen_chain, cutoff=cutoff)

            for res_id, res_name in ag_bridges:
                key = f"{antigen_chain}:{res_id + uniprot_correction}-{res_name}"
                bridge_counter[key] = bridge_counter.get(key, 0) + 1
                residue_counter[(antigen_chain, res_id, res_name)] += 1

            for (res_id, res_name), partner in pairs:
                ag_key = (antigen_chain, res_id, res_name)
                partner_counter.setdefault(ag_key, Counter())[partner] += 1

            hc_seen.update(hc_bridges)
            lc_seen.update(lc_bridges)

        for key in hc_seen:
            hc_bridge_counter[key] = hc_bridge_counter.get(key, 0) + 1

        for key in lc_seen:
            lc_bridge_counter[key] = lc_bridge_counter.get(key, 0) + 1

    def antibody_lines(counter, chain):
        """Frequency lines for one antibody chain, tagged with the CDR each residue sits in."""
        out = [f"\n{chain} bridge residue frequencies:"]
        for (res_id, res_name), count in sorted(counter.items()):
            cdr = get_cdr_name(res_id, cdr_ranges, chain)
            region = cdr if cdr else "framework"
            out.append(f"  {res_id}-{res_name} [{region}]: {count}")
        return out

    lines = [f"\n{'='*60}", f"Section: {label} | cutoff={cutoff}Å | {n} structures", f"{'='*60}"]
    lines.append("\nAntigen bridge residue frequencies:")
    for key, count in sorted(bridge_counter.items()):
        lines.append(f"  {key}: {count}")

    if hc_bridge_counter:
        lines += antibody_lines(hc_bridge_counter, 'HC')

    if lc_bridge_counter:
        lines += antibody_lines(lc_bridge_counter, 'LC')

    for line in lines:
        print(line)

    return lines, residue_counter, partner_counter


def count_kept_bridges(subset, kept_partners, cutoff=4.0):
    """Number of structures in which each antigen residue bridges with at least one kept partner."""
    counter = Counter()
    for sample in subset['sample_name']:
        cif = os.path.join(CIF_DIR, f"{sample}.cif")
        for antigen_chain, _ in ANTIGEN_CHAINS:
            *_, pairs = find_salt_bridges(cif, antigen_chain, cutoff=cutoff)
            # A set, so a structure with two kept partners still counts once.
            bridged = {(antigen_chain, res_id, res_name)
                       for (res_id, res_name), partner in pairs
                       if partner in kept_partners.get((antigen_chain, res_id, res_name), {})}
            counter.update(bridged)
    return counter


def salt_bridge_table(section_counts, partners, top_n):
    """
    One row per antigen residue that bridges at 4 Å in any prediction, with its
    count in every section and its kept antibody partners at 4 Å over all predictions.
    A residue without a kept partner has 0 counts and no partners.
    """
    corrections = dict(ANTIGEN_CHAINS)
    residues = sorted(partners)

    rows = []
    for chain, res_id, res_name in residues:
        key = (chain, res_id, res_name)
        row = {
            "chain": chain,
            "residue_id": res_name,
            "residue_num": res_id,
            "residue_num_uniprot": res_id + corrections[chain],
        }
        for subset, suffix in (("all", "all"), ("top", f"top{top_n}_{RANK_COLUMN}")):
            for cutoff in CUTOFFS:
                row[f"salt_bridges_{cutoff:g}a_{suffix}"] = section_counts[(subset, cutoff)].get(key, 0)
        # Partners as 'D59(H):514; D1(L):12', most frequent first.
        ranked = sorted(partners[key].items(), key=lambda item: (-item[1], item[0]))
        row[f"partners_{CUTOFFS[0]:g}a_all"] = "; ".join(
            f"{AA3_TO_1.get(p_name, p_name)}{p_id}({p_chain}):{count}"
            for (p_chain, p_id, p_name), count in ranked)
        rows.append(row)

    return pd.DataFrame(rows)


if __name__ == "__main__":
    # ── 1. Arguments and antigen chains ─────────────
    parser = argparse.ArgumentParser(description="Visualize salt bridge frequencies across top structures")
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
                        help="Number of top predictions by iptm for the top-N sections")
    parser.add_argument("--partner_fraction", type=float, default=0.10,
                        help="Min fraction of all predictions in which an antibody partner "
                             "bridges at 4 Å to be kept in the table")
    args = parser.parse_args()

    ANTIGEN_CHAINS = [(args.antigen_chain, args.uniprot_correction)]   ### [(chain_id, uniprot_correction), ...]
    if args.antigen_chain_2:
        ANTIGEN_CHAINS.append((args.antigen_chain_2, args.uniprot_correction_2))
    print(f"Antigen chains: {ANTIGEN_CHAINS}")

    # ── 2. Paths ─────────────
    structure_name = str(args.structure_name).rstrip("/")
    name = os.path.basename(structure_name)
    sorted_output = f'{structure_name}/full_antigen/sorted_output'
    csv_path      = f'{sorted_output}/csv_files_metrics/{name}.csv'
    analysis_dir  = f'{structure_name}/full_antigen/analysis_results'
    CIF_DIR       = f'{sorted_output}/cif_folder/{name}_cif'

    if not os.path.exists(csv_path):
        raise SystemExit(f"ERROR: output CSV not found: {csv_path}")

    print(f"\n=== Processing structure: {name} ===")
    df = pd.read_csv(csv_path)

    # ── 3. CDR ranges ─────────────
    # CDR ranges are derived once as the antibody is the same.
    first_cif = os.path.join(CIF_DIR, f"{df['sample_name'].iloc[0]}.cif")
    cdr_ranges = find_cdr_ranges(first_cif)

    all_lines = format_cdr_ranges(cdr_ranges)
    for line in all_lines:
        print(line)

    # ── 4. Salt bridges per 'subset x cutoff' ─────────────
    subsets = {
        "all": (df.reset_index(drop=True), f"all_{len(df)}"),
        "top": (df.nlargest(args.top_n, RANK_COLUMN).reset_index(drop=True),
                f"top{args.top_n}_{RANK_COLUMN}"),
    }

    section_counts = {}
    partners = {}
    for subset, (subset_df, label) in subsets.items():
        for cutoff in CUTOFFS:
            lines, residue_counter, partner_counter = count_salt_bridges(
                subset_df, cdr_ranges, cutoff=cutoff, label=label)
            all_lines += lines
            section_counts[(subset, cutoff)] = residue_counter
            if subset == "all" and cutoff == CUTOFFS[0]:
                partners = partner_counter

    # Keep the partners bridging at 4 Å in more than 'partner_fraction' of all predictions,
    # then recount every section from those partners alone. The .txt report stays unfiltered.
    kept_partners = {
        key: {partner: count for partner, count in counts.items()
              if count > args.partner_fraction * len(df)}
        for key, counts in partners.items()
    }
    kept_counts = {
        (subset, cutoff): count_kept_bridges(subset_df, kept_partners, cutoff=cutoff)
        for subset, (subset_df, _) in subsets.items()
        for cutoff in CUTOFFS
    }

    # ── 5. Write the report and the table ─────────────
    os.makedirs(analysis_dir, exist_ok=True)
    out_path = os.path.join(analysis_dir, f"salt_bridge_frequencies_{name}.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(all_lines) + "\n")
    print(f"Salt bridge frequencies written to: {out_path}")

    table = salt_bridge_table(kept_counts, kept_partners, args.top_n)
    table_path = os.path.join(analysis_dir, f"salt_bridge_frequencies_{name}.csv")
    table.to_csv(table_path, index=False)
    print(f"Salt bridge table written to: {table_path}")

