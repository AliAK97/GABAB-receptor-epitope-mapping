"""
Extract structural and AlphaFold3 confidence metrics for antibody-antigen predictions.

For every predicted seed-sample this script collects:

  - distance metrics (RMSD, GDT_TS, lDDT) against a reference structure.
  - the predicted epitope residues at 4, 3 and 2.5 A contact cutoffs;
  - the AlphaFold3 confidence scores (ptm, iptm, ranking score).

Results are merged into one CSV per complex.

One or two antigen chains are supported. The second chain is optional and is
enabled by an 'antigen_chain_2' column in the metadata CSV; all per-antigen
outputs are then duplicated with an '_ag2' suffix.
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder, NeighborSearch, Superimposer
from Bio.Align import PairwiseAligner
from functools import reduce


CONTACT_CUTOFFS = [(4.0, "4a"), (3.0, "3a"), (2.5, "2.5a")]   ### (cutoff in Å, column tag) for the epitope


def extract_chain_sequence(structure, chain_id):
    """One-letter sequence of a chain, joined over its peptide fragments."""
    chain = structure[0][chain_id]
    ppb = PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(chain):
        sequence += str(pp.get_sequence())
    return sequence


def align_sequences(seq1, seq2):
    """Best global pairwise alignment of two sequences."""
    aligner = PairwiseAligner()
    aligner.match = 5
    aligner.mismatch = 0
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5
    return aligner.align(seq1, seq2)[0]


def get_matching_atoms(aln, model_chain, native_chain):
    """CA atoms of the residues paired by the alignment, as (model atoms, native atoms)."""
    model_residues = list(model_chain.get_residues())
    native_residues = list(native_chain.get_residues())

    model_atoms = []
    native_atoms = []

    model_aligned = aln.aligned[0]
    native_aligned = aln.aligned[1]

    for (m_start, m_end), (n_start, n_end) in zip(model_aligned, native_aligned):
        for mi, ni in zip(range(m_start, m_end), range(n_start, n_end)):
            if mi < len(model_residues) and ni < len(native_residues):
                m_res = model_residues[mi]
                n_res = native_residues[ni]
                if 'CA' in m_res and 'CA' in n_res:
                    model_atoms.append(m_res['CA'])
                    native_atoms.append(n_res['CA'])

    return model_atoms, native_atoms


def get_coords_from_atoms(atom_list):
    """Coordinates of a list of atoms as an (n, 3) array."""
    return np.array([atom.get_coord() for atom in atom_list])


def calc_rmsd(coords1, coords2):
    """RMSD between two matched coordinate sets, without superposition."""
    diff = coords1 - coords2
    return np.sqrt((diff * diff).sum() / len(coords1))


def calc_gdt_ts(coords1, coords2):
    """GDT-TS: mean percentage of atoms within 1, 2, 4 and 8 Å."""
    if len(coords1) != len(coords2):
        print("Unequal length of residues")
        return None
    distances = np.sqrt(np.sum((coords1 - coords2) ** 2, axis=1))
    p1 = np.mean(distances <= 1.0) * 100
    p2 = np.mean(distances <= 2.0) * 100
    p4 = np.mean(distances <= 4.0) * 100
    p8 = np.mean(distances <= 8.0) * 100
    return (p1 + p2 + p4 + p8) / 4.0


def calc_lddt(coords1, coords2, cutoff=15.0,
              threshold_values=None):
    """
    lDDT: percentage of atom pairs within 'cutoff' Å in the reference whose distance
    is kept within 0.5, 1, 2 and 4 Å in the model, averaged over the four thresholds.
    """
    if threshold_values is None:
        threshold_values = [0.5, 1.0, 2.0, 4.0]
    if len(coords1) != len(coords2):
        print("Unequal length of residues")
        return None

    n_atoms = len(coords1)

    dist_ref = np.linalg.norm(coords1[:, None, :] - coords1[None, :, :], axis=-1)
    dist_pred = np.linalg.norm(coords2[:, None, :] - coords2[None, :, :], axis=-1)

    preserved = []
    for t in threshold_values:
        mask = (dist_ref < cutoff)
        diff = np.abs(dist_ref - dist_pred)
        preserved_pairs = np.sum((diff <= t) & mask) - n_atoms
        total_pairs = np.sum(mask) - n_atoms
        preserved.append(preserved_pairs / total_pairs if total_pairs > 0 else 0.0)

    return np.mean(preserved) * 100


def parse_af3_json_metrics(json_path, antigen_id, heavy_id, light_id, antigen2_id=None):
    """AF3 confidence scores of every summary .json file, one row per sample."""
    full_results = {}
    for file in os.listdir(json_path):
        if not file.endswith(".json"):
            continue

        full_file_path = os.path.join(json_path, file)
        sample_name = os.path.splitext(file)[0]
        with open(full_file_path, 'r') as f:
            data = json.load(f)

        chain_ids = [antigen_id, heavy_id, light_id]
        if antigen2_id:
            chain_ids.append(antigen2_id)
        n_chains = len(chain_ids)
        sorted_chain_ids = sorted(chain_ids)
        chain_index_map = {chain: idx for idx, chain in enumerate(sorted_chain_ids)}

        role_map = {antigen_id: "antigen", heavy_id: "heavy", light_id: "light"}
        if antigen2_id:
            role_map[antigen2_id] = "antigen2"
        index_to_role = {chain_index_map[c]: role_map[c] for c in chain_ids}

        result = {}

        result["sample_name"] = sample_name

        for i, val in enumerate(data["chain_ptm"]):
            role = index_to_role[i]
            result[f"{role}_ptm"] = val

        for i, val in enumerate(data["chain_iptm"]):
            role = index_to_role[i]
            result[f"{role}_iptm"] = val

        for i in range(n_chains):
            for j in range(n_chains):
                role_i = index_to_role[i]
                role_j = index_to_role[j]
                val = data["chain_pair_iptm"][i][j]

                if i == j:
                    result[f"{role_i}_pair_iptm"] = val
                else:
                    pair_name = "_".join(sorted([role_i, role_j]))
                    result[f"{pair_name}_pair_iptm"] = val

        for i in range(n_chains):
            for j in range(n_chains):
                role_i = index_to_role[i]
                role_j = index_to_role[j]
                val = data["chain_pair_pae_min"][i][j]

                if i == j:
                    result[f"{role_i}_pair_pae_min"] = val
                else:
                    pair_name = "_".join(sorted([role_i, role_j]))
                    result[f"{pair_name}_pair_pae_min"] = val

        result["iptm"] = data["iptm"]
        result["ptm"] = data["ptm"]
        result["ranking_score"] = data["ranking_score"]
        result["mean_iptm"] = np.mean(data["chain_iptm"])
        result["mean_ptm"] = np.mean(data["chain_ptm"])

        full_results[sample_name] = result
    return pd.DataFrame.from_dict(full_results, orient="index").reset_index(drop=True)


def distance_metrics_dict_to_dataframe(results, dataset_name):
    """Distance metrics as one row per sample, with columns <dataset>_<metric>_<chain>."""
    metrics = ["rmsd", "lddt", "gdt_ts"]
    chain_types = ["overall", "heavy_chain", "light_chain", "antigen_chain", "antigen_chain2"]
    records = []

    for sample, chains in results.items():
        for chain_type in chain_types:
            for metric in metrics:
                value = chains.get(chain_type, {}).get(metric)
                if value is not None:
                    records.append({
                        "sample_name": sample,
                        "Chain": chain_type,
                        "Metric": metric,
                        "Value": value
                    })

    df = pd.DataFrame(records)
    df = df.sort_values(by=["sample_name", "Metric", "Chain"])
    df = df.reset_index(drop=True)

    rows = []

    for _, group in df.groupby("sample_name"):
        group = group.copy()
        group["column_name"] = dataset_name + "_" + group["Metric"] + "_" + group["Chain"]

        wide = group.pivot(index="sample_name", columns="column_name", values="Value")
        wide.columns.name = None
        wide = wide.reset_index()
        rows.append(wide)

    final_df = pd.concat(rows, ignore_index=True)
    return final_df


class PredictionEvaluator:
    """Loads structures and compares a predicted structure to the experimental one."""
    def __init__(self):
        self.cif_parser = MMCIFParser(QUIET=True)
        self.pdb_parser = PDBParser(QUIET=True)

    def read_csv_file(self, structure_name, csv_path):
        """Heavy, light and antigen chain IDs of the structure from the metadata CSV."""
        df = pd.read_csv(csv_path)
        if structure_name not in df['structure_name'].values:
            raise ValueError(
                f"No row for structure '{structure_name}' in the metadata file")

        heavy_chain_id = df.loc[df['structure_name'] == structure_name, 'heavy_chain'].values[0]
        light_chain_id = df.loc[df['structure_name'] == structure_name, 'light_chain'].values[0]
        antigen_id = df.loc[df['structure_name'] == structure_name, 'antigen_chain'].values[0]

        antigen2_id = None
        if 'antigen_chain_2' in df.columns:
            val = df.loc[df['structure_name'] == structure_name, 'antigen_chain_2'].values[0]
            if pd.notna(val) and val != '':
                antigen2_id = val

        return heavy_chain_id, light_chain_id, antigen_id, antigen2_id

    def load_structure(self, file_path):
        """Structure from a .pdb or .cif file."""
        if file_path.endswith(".pdb"):
            structure = self.pdb_parser.get_structure("model", file_path)
        else:
            structure = self.cif_parser.get_structure("model", file_path)
        return structure

    def calculate_metrics(self, native_structure, model_structure,
                          target_chains, heavy_chain_id=None,
                          light_chain_id=None, ag_chain_id=None,
                          ag_chain2_id=None, verbose=False):
        """
        RMSD, GDT-TS and lDDT of the model against the native structure, overall and
        per chain, after superimposing the 'target_chains' (the antigen) by CA atoms.
        """
        model_chains = {chain.id for chain in model_structure[0]}
        native_chains = {chain.id for chain in native_structure[0]}
        common_chains = model_chains & native_chains

        target_chains = set(target_chains)

        results = {
            "overall": {},
            "heavy_chain": {},
            "light_chain": {},
            "antigen_chain": {},
            "antigen_chain2": {}
        }

        all_model_atoms = []
        all_native_atoms = []
        super_model_atoms = []
        super_native_atoms = []
        per_chain_data = {}

        for chain_id in common_chains:
            model_seq = extract_chain_sequence(model_structure, chain_id)
            native_seq = extract_chain_sequence(native_structure, chain_id)

            aln = align_sequences(model_seq, native_seq)
            model_chain = model_structure[0][chain_id]
            native_chain = native_structure[0][chain_id]

            model_atoms, native_atoms = get_matching_atoms(aln, model_chain, native_chain)

            if not model_atoms:
                continue

            all_model_atoms.extend(model_atoms)
            all_native_atoms.extend(native_atoms)

            if chain_id in target_chains:
                super_model_atoms.extend(model_atoms)
                super_native_atoms.extend(native_atoms)

            per_chain_data[chain_id] = (model_atoms, native_atoms)

        sup = Superimposer()
        sup.set_atoms(super_native_atoms, super_model_atoms)
        sup.apply(model_structure.get_atoms())

        coords_model_all = get_coords_from_atoms(all_model_atoms)
        coords_native_all = get_coords_from_atoms(all_native_atoms)

        results["overall"]["rmsd"] = calc_rmsd(coords_native_all, coords_model_all)
        results["overall"]["gdt_ts"] = calc_gdt_ts(coords_native_all, coords_model_all)
        results["overall"]["lddt"] = calc_lddt(coords_native_all, coords_model_all)

        chain_labels = [
            ("heavy_chain", heavy_chain_id),
            ("light_chain", light_chain_id),
            ("antigen_chain", ag_chain_id),
            ("antigen_chain2", ag_chain2_id),
        ]
        for label, cid in chain_labels:
            if cid and cid in per_chain_data:
                model_atoms, native_atoms = per_chain_data[cid]
                coords_model = get_coords_from_atoms(model_atoms)
                coords_native = get_coords_from_atoms(native_atoms)

                results[label]["rmsd"] = calc_rmsd(coords_native, coords_model)
                results[label]["gdt_ts"] = calc_gdt_ts(coords_native, coords_model)
                results[label]["lddt"] = calc_lddt(coords_native, coords_model)

        if verbose:
            print("\nOverall Structure Metrics:")
            print(f"RMSD: {results['overall']['rmsd']:.3f} Å")
            print(f"GDT_TS: {results['overall']['gdt_ts']:.2f}%")
            print(f"LDDT: {results['overall']['lddt']:.2f}%")

            for label, cid in chain_labels:
                if cid and "rmsd" in results[label]:
                    print(f"\n{label.replace('_', ' ').title()} ({cid}) Metrics:")
                    print(f"RMSD: {results[label]['rmsd']:.3f} Å")
                    print(f"GDT_TS: {results[label]['gdt_ts']:.2f}%")
                    print(f"LDDT: {results[label]['lddt']:.2f}%")
        return results


class EpitopeFinder:
    """Contacts between the antigen and antibody chains of one predicted structure."""
    def __init__(self, structure_path, antigen_chain, antibody_chains, distance_cutoff=4.0):
        self.structure_path = structure_path
        self.antigen_chain = antigen_chain
        self.antibody_chains = antibody_chains
        self.distance_cutoff = distance_cutoff

        if structure_path.lower().endswith(".cif"):
            parser = MMCIFParser(QUIET=True, auth_residues=False)
        elif structure_path.lower().endswith(".pdb"):
            parser = PDBParser(QUIET=True)
        else:
            raise ValueError(f"Unsupported file format: {structure_path}")

        self.structure = parser.get_structure("complex", structure_path)

        # Atoms of the antigen and antibody chains in the first model.
        model = self.structure[0]
        self.antigen_atoms = [atom for atom in model[self.antigen_chain].get_atoms()]
        self.antibody_atoms = [
            atom for ch in self.antibody_chains for atom in model[ch].get_atoms()
        ]

        # Spatial index of the antibody atoms.
        self.ns = NeighborSearch(self.antibody_atoms)

    def compute_contacts(self):
        """For each antigen residue, the number of antibody residues within 'distance_cutoff' Å."""
        results = []
        for residue in self.structure[0][self.antigen_chain]:
            if not residue.has_id("CA"):
                continue
            
            contact_residues = set()
            for atom in residue.get_atoms():
                for close_atom in self.ns.search(atom.coord, self.distance_cutoff):
                    contact_residues.add(close_atom.get_parent())

            contact_count = len(contact_residues)

            results.append({
                "res_id": residue.get_id()[1],   ### CIF residue number, starting at 1
                "res_name": residue.get_resname(),
                "contact_count": contact_count,
                "classification": self._classify(contact_count)
            })

        return results

    def _classify(self, contact_count):
        """'contact' if the residue touches any antibody residue, else 'noncontact'."""
        if contact_count == 0:
            return "noncontact"
        return "contact"

    def epitope_residue_list(self, results):
        """Residue IDs and names of the antigen residues in contact."""
        return [
            {"res_id": int(r["res_id"]), "res_name": r["res_name"]}
            for r in results
            if r["classification"] == "contact"
        ]

    def compute_residue_min_distances(self):
        """
        For each antigen residue, the minimum atom-to-atom distance to any antibody
        atom. Used to draw the epitope heatmap on the antigen chain.
        """
        ab_coords = np.array([a.get_coord() for a in self.antibody_atoms])
        residues = [r for r in self.structure[0][self.antigen_chain].get_residues()
                    if r.has_id("CA")]
        min_dists = []
        for res in residues:
            ag_coords = np.array([a.get_coord() for a in res.get_atoms()])
            diffs = ag_coords[:, None, :] - ab_coords[None, :, :]
            min_dist = np.sqrt((diffs ** 2).sum(axis=-1)).min()
            min_dists.append(float(min_dist))
        return min_dists


def analyse_epitope(cif_path, antigen_chain, antibody_chains, cutoff):
    """Contact analysis of one antigen chain at one distance cutoff."""
    finder = EpitopeFinder(
        cif_path,
        antigen_chain=antigen_chain,
        antibody_chains=antibody_chains,
        distance_cutoff=cutoff)
    results = finder.compute_contacts()

    residues = finder.epitope_residue_list(results)
    contact_count = sum(
        r["contact_count"] for r in results
        if r["classification"] == "contact"
    )

    return finder, {
        "residues": residues,
        "length": len(residues),
        "contact_count": contact_count,
    }


def epitope_columns(cif_path, antibody_chains, antigen_specs):
    """Epitope and minimum-distance columns of one predicted sample."""
    columns = {}
    finders = {}

    for suffix, chain_id in antigen_specs:
        for cutoff, tag in CONTACT_CUTOFFS:
            finder, info = analyse_epitope(cif_path, chain_id, antibody_chains, cutoff)
            finders.setdefault(suffix, finder)

            # 4 Å keeps the bare column names; the other cutoffs are tagged.
            cutoff_tag = "" if tag == "4a" else f"_{tag}"
            columns[f"contact_residues{suffix}{cutoff_tag}"] = info["residues"]
            columns[f"length_epitope_predicted_{tag}{suffix}"] = info["length"]
            columns[f"contact_count_pred_{tag}{suffix}"] = info["contact_count"]

    for suffix, chain_id in antigen_specs:
        prefix = "ag_res" if suffix == "" else "ag2_res"
        min_dists = finders[suffix].compute_residue_min_distances()
        for i, distance in enumerate(min_dists):
            columns[f"{prefix}_{i + 1}_min_dist"] = distance

    return columns


def run_function(structure_name, exp_structure_file):
    """
    Metrics of every seed-sample of one structure, compared to the experimental
    structure 'exp_structure_file' (e.g. 6UO8_chainB.cif), written to one CSV file.
    """
    evaluator = PredictionEvaluator()
    ag_aligned_metrics = {}
    epitope_results = []

    # ── 1. Paths ─────────────
    structure_name = str(structure_name).rstrip("/")
    name = os.path.basename(structure_name)
    sorted_output = f'{structure_name}/full_antigen/sorted_output'
    meta_csv_path = f'{sorted_output}/meta_data.csv'
    exp_structure_path = f'{sorted_output}/experimental_structure/{exp_structure_file}'
    pred_cif_folder  = f'{sorted_output}/cif_folder/{name}_cif'
    pred_json_folder = f'{sorted_output}/json_folder/{name}_json_summary'

    # ── 2. Chain IDs and the experimental structure ─────────────
    heavy_chain_id, light_chain_id, antigen_id, antigen2_id = evaluator.read_csv_file(name, meta_csv_path)
    print(f'heavy: {heavy_chain_id}, light: {light_chain_id}, antigen: {antigen_id}, antigen2: {antigen2_id}')
    exp_structure = evaluator.load_structure(exp_structure_path)

    antigen_specs = [("", antigen_id)]
    if antigen2_id:
        antigen_specs.append(("_ag2", antigen2_id))

    # ── 3. Metrics of every predicted structure ─────────────
    for file in os.listdir(pred_cif_folder):
        if file.endswith(".cif"):
            full_file_path = os.path.join(pred_cif_folder, file)
            sample_name = os.path.splitext(file)[0]
            print(f" ++++++++++++ Processing file: {sample_name}")
            pred_structure = evaluator.load_structure(full_file_path)

            # RMSD, GDT-TS and lDDT after aligning the antigen to the experimental structure.
            print(" ------------ Metrics for antigen aligned structure")
            target_chains = [antigen_id] + ([antigen2_id] if antigen2_id else [])
            ag_aligned_metrics_single = evaluator.calculate_metrics(
                exp_structure,
                pred_structure,
                target_chains=target_chains,
                heavy_chain_id=heavy_chain_id,
                light_chain_id=light_chain_id,
                ag_chain_id=antigen_id,
                ag_chain2_id=antigen2_id,
                verbose=True)

            ag_aligned_metrics[sample_name] = ag_aligned_metrics_single

            epitope_results.append({
                "sample_name": sample_name,
                **epitope_columns(full_file_path,
                                  antibody_chains=(heavy_chain_id, light_chain_id),
                                  antigen_specs=antigen_specs),
            })

    # ── 4. Merge and write the CSV ─────────────
    df_af3_metrics = parse_af3_json_metrics(pred_json_folder,
                                            antigen_id,
                                            heavy_chain_id,
                                            light_chain_id,
                                            antigen2_id)
    df_epitope = pd.DataFrame(epitope_results)
    ag_align_df = distance_metrics_dict_to_dataframe(ag_aligned_metrics, "ag_align")

    dfs_list = [df_af3_metrics, df_epitope, ag_align_df]
    merged_df = reduce(lambda left, right: pd.merge(left, right, on="sample_name"), dfs_list)

    merged_df.to_csv(f'{sorted_output}/csv_files_metrics/{name}.csv', index=False)
    return merged_df


def main(structure_name, exp_structure_file):
    """Runs the extraction unless the structure is not sorted yet or its CSV already exists."""
    structure_name = str(structure_name).rstrip("/")
    name = os.path.basename(structure_name)
    sorted_output = f'{structure_name}/full_antigen/sorted_output'
    folder_csvs = f'{sorted_output}/csv_files_metrics'
    folder_cifs = f'{sorted_output}/cif_folder'
    os.makedirs(folder_csvs, exist_ok=True)

    csv_ids ={os.path.splitext(f)[0] for f in os.listdir(folder_csvs) if f.endswith(".csv")}
    sorted_ids = {f.replace("_cif", "") for f in os.listdir(folder_cifs)
                  if not f.startswith('.')}

    if name not in sorted_ids:
        print(f"{name} is not yet sorted (no {name}_cif in {folder_cifs}). Skipping.")
        return

    if name in csv_ids:
        print(f"Skipping {name}, CSV already exists.")
        return

    print(f"Processing {name}")
    try:
        run_function(structure_name, exp_structure_file)
    except Exception as e:
        print(f"Error processing {name}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract metrics from AlphaFold3 antibody-antigen predictions."
    )
    parser.add_argument("structure_name",
                        help="Path to the structure directory")
    parser.add_argument("--exp_structure", type=str, required=True,
                        help="File name of the experimental structure in "
                             "sorted_output/experimental_structure/ (e.g. 6UO8_chainB.cif)")
    args = parser.parse_args()

    main(args.structure_name, args.exp_structure)
