"""
Computational alanine scan of the antigen-antibody interface with PyRosetta.
Each of the top N ipTM-ranked complexes is relaxed with Cartesian FastRelax, then
every antigen residue within the interface cutoff is mutated to alanine. Its binding
ΔΔG (REU) is the mean mutant minus the mean wild-type binding energy over the
packing trials. Ala, Gly, Pro and disulfide-bonded Cys are not scanned.

Writes into analysis_results/:
  alanine_scan_<name>_<chain>_chunk<N>.csv
"""


import argparse
import glob
import os
import pandas as pd
import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.protocols.scoring import Interface


ANTIBODY_CHAINS = "HL"
SKIP_RESIDUES = {"A", "G", "P"}  ### Ala is already Ala; Gly and Pro to Ala change the backbone
RANK_COLUMN = "iptm"


def relax_pose(pose, scorefxn, max_iter, nstruct=10):
    """Cartesian FastRelax replicates of the pose and returns the lowest-energy one."""
    replicates = []
    for _ in range(nstruct):
        trial = pose.clone()
        fr = FastRelax()
        fr.set_scorefxn(scorefxn)
        fr.max_iter(max_iter)
        fr.cartesian(True)
        fr.min_type("lbfgs_armijo_nonmonotone")
        fr.apply(trial)
        replicates.append((scorefxn(trial), trial))

    return min(replicates, key=lambda replicate: replicate[0])[1]


def get_interface_residues(pose, antigen_chain, ref_scorefxn, cutoff=12.0):
    """Pose indices of the antigen residues within 'cutoff' Å of the antibody."""
    ref_scorefxn(pose)

    interface = Interface(1)
    interface.distance(cutoff)
    interface.calculate(pose)

    interface_indices = []
    for rosetta_idx in range(1, pose.total_residue() + 1):
        if pose.pdb_info().chain(rosetta_idx) != antigen_chain:
            continue
        if interface.is_interface(rosetta_idx):
            interface_indices.append(rosetta_idx)

    return interface_indices


def calc_binding_energy(pose, scorefxn, partners):
    """Binding energy (dG_separated) of the complex, repacked before and after separation."""
    test_pose = pose.clone()

    ia = InterfaceAnalyzerMover(partners)
    ia.set_scorefunction(scorefxn)
    ia.set_pack_input(True)
    ia.set_pack_separated(True)
    ia.apply(test_pose)

    return ia.get_interface_dG()


def mutate_to_ala(pose, rosetta_idx):
    """Copy of the pose with one residue mutated to alanine."""
    mut_pose = pose.clone()
    result = mutate_residue(mut_pose, rosetta_idx, 'A')
    if result is not None:
        mut_pose = result

    if mut_pose.residue(rosetta_idx).name1() != 'A':
        raise RuntimeError(
            f"mutate_residue failed to mutate residue {rosetta_idx} to Ala "
            f"(still {mut_pose.residue(rosetta_idx).name1()})."
        )
    return mut_pose


def compute_ddg(pose, rosetta_idx, scorefxn, partners, trials=5):
    """Mean mutant minus mean wild-type binding energy over 'trials' packing trials."""
    wt_base = pose.clone()
    mut_base = mutate_to_ala(pose, rosetta_idx)

    wt_energies = []
    mut_energies = []
    for _ in range(trials):
        wt_energies.append(
            calc_binding_energy(wt_base.clone(), scorefxn, partners)
        )
        mut_energies.append(
            calc_binding_energy(mut_base.clone(), scorefxn, partners)
        )

    mean_wt = sum(wt_energies) / len(wt_energies)
    mean_mut = sum(mut_energies) / len(mut_energies)
    return mean_mut - mean_wt


def run_alanine_scan(pdb_path, antigen_chain, partners, uniprot_correction,
                     cutoff=12.0, trials=5, relax_iter=5, relax_nstruct=10,
                     seed=42):
    """ΔΔG of every scanned antigen interface residue in one structure."""
    if seed != 0:
        pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
        print(f"Rosetta RNG seed set to: {seed}")
    else:
        print("WARNING: No fixed seed. Results will differ between runs.")

    ddg_sf = pyrosetta.create_score_function("ref2015_cart")   ### for the Cartesian relax
    ref_sf = pyrosetta.create_score_function("ref2015")        ### for the binding energies

    pose = pose_from_pdb(pdb_path)

    movable_jumps = pyrosetta.rosetta.utility.vector1_int()
    movable_jumps.append(1)
    setup_foldtree(pose, partners, movable_jumps)

    print(f"  Relaxing {os.path.basename(pdb_path)} "
          f"({relax_nstruct} cartesian replicates)...")
    pose = relax_pose(pose, ddg_sf, max_iter=relax_iter, nstruct=relax_nstruct)

    results = []

    interface_rosetta_indices = get_interface_residues(
        pose, antigen_chain, ref_sf, cutoff=cutoff
    )

    for rosetta_idx in interface_rosetta_indices:
        pdb_res_num = pose.pdb_info().number(rosetta_idx)
        wt_aa = pose.residue(rosetta_idx).name1()
        uniprot_res_num = pdb_res_num + uniprot_correction
        mutation_label = f"{wt_aa}{uniprot_res_num}A"

        if wt_aa in SKIP_RESIDUES:
            print(f"  {wt_aa}{pdb_res_num}A: skipped (Ala, Gly or Pro)")
            continue

        res = pose.residue(rosetta_idx)
        if res.name1() == 'C' and res.has_variant_type(
                pyrosetta.rosetta.core.chemical.VariantType.DISULFIDE):
            print(f"  {wt_aa}{pdb_res_num}A: skipped (disulfide-bonded)")
            continue

        mean_ddg = compute_ddg(pose, rosetta_idx, ref_sf, partners, trials)

        results.append({
            "mutation": mutation_label,
            "ddG_REU":  round(mean_ddg, 3),
        })
        print(f"  {wt_aa}{pdb_res_num}A: ddG_bind = {mean_ddg:.3f} REU")

    return pd.DataFrame(results)


if __name__ == "__main__":
    pyrosetta.init("-mute all")

    # ── 1. Arguments ─────────────
    parser = argparse.ArgumentParser(
        description="Rigid-backbone computational alanine scanning (Kortemme & Baker 2002)"
    )
    parser.add_argument("structure_name",
                        help="Path to the structure directory")
    parser.add_argument("--antigen_chain",       type=str, required=True,
                        help="Chain ID of the antigen (e.g. C)")
    parser.add_argument("--antigen_chain_2",     type=str, default="",
                        help="Chain ID of a second antigen chain. Scanned as its own "
                             "run against the antibody, not merged with the first.")
    parser.add_argument("--uniprot_correction",  type=int, default=0,
                        help="Offset added to PDB residue numbers for UniProt numbering")
    parser.add_argument("--uniprot_correction_2", type=int, default=0,
                        help="Offset for --antigen_chain_2")
    parser.add_argument("--top_n",               type=int, default=30,
                        help="Number of top predicted structures to scan (default 30)")
    parser.add_argument("--trials",              type=int, default=5,
                        help="Packing trials per residue per structure (default 5, min 3)")
    parser.add_argument("--cutoff",              type=float, default=12.0,
                        help="Interface distance cutoff in Angstrom (default 12.0)")
    parser.add_argument("--seed",                type=int, default=42,
                        help="Rosetta RNG seed (default 42; 0 = disabled)")
    parser.add_argument("--relax_iter",          type=int, default=5,
                        help="FastRelax max_iter (default 5)")
    parser.add_argument("--relax_nstruct",       type=int, default=10,
                        help="Number of FastRelax replicates; best scoring kept (default 10)")
    parser.add_argument("--chunk_index",         type=int, default=-1,
                        help="0-based chunk index for SLURM array jobs (-1 = run all)")
    parser.add_argument("--chunk_size",          type=int, default=3,
                        help="Number of structures per chunk (default 3)")
    args = parser.parse_args()

    # ── 2. Paths and the structures of this chunk ─────────────
    structure_name = str(args.structure_name).rstrip("/")
    name = os.path.basename(structure_name)
    sorted_output = f'{structure_name}/full_antigen/sorted_output'
    csv_path      = f'{sorted_output}/csv_files_metrics/{name}.csv'
    pdb_folder    = f'{sorted_output}/pdb_folder'
    analysis_dir  = f'{structure_name}/full_antigen/analysis_results'
    os.makedirs(analysis_dir, exist_ok=True)

    df_scores = pd.read_csv(csv_path)
    list_best = df_scores.nlargest(args.top_n, RANK_COLUMN)["sample_name"].tolist()

    if args.chunk_index >= 0:
        start     = args.chunk_index * args.chunk_size
        end       = start + args.chunk_size
        list_best = list_best[start:end]
        print(f"Chunk {args.chunk_index}: structures {start}–{end - 1}: {list_best}")
    else:
        print(f"Top {args.top_n} structures selected: {list_best}")

    pdb_by_stem = {os.path.splitext(os.path.basename(path))[0]: path
                   for path in glob.glob(os.path.join(pdb_folder, "*", "*.pdb"))}

    # partners "A_HL", then "B_HL"
    antigen_chains = [(args.antigen_chain, args.uniprot_correction)]
    if args.antigen_chain_2:
        antigen_chains.append((args.antigen_chain_2, args.uniprot_correction_2))

    chunk_tag = f"_chunk{args.chunk_index}" if args.chunk_index >= 0 else ""

    # ── 3. Alanine scan per antigen chain ─────────────
    summary = []
    for antigen_chain, uniprot_correction in antigen_chains:
        partners = f"{antigen_chain}_{ANTIBODY_CHAINS}"
        print(f"\n{'#'*60}")
        print(f"# Antigen chain {antigen_chain}  (partners: {partners})")
        print(f"{'#'*60}")

        all_results = []
        missing, failed = [], []

        for sample in list_best:
            pdb_path = pdb_by_stem.get(sample)
            if pdb_path is None:
                print(f"\n  ERROR: no PDB file found for {sample}")
                missing.append(sample)
                continue

            print(f"\n{'='*60}")
            print(f"Processing: {os.path.basename(pdb_path)}")
            print(f"{'='*60}")

            try:
                scan_df = run_alanine_scan(
                    pdb_path=pdb_path,
                    antigen_chain=antigen_chain,
                    partners=partners,
                    uniprot_correction=uniprot_correction,
                    cutoff=args.cutoff,
                    trials=args.trials,
                    relax_iter=args.relax_iter,
                    relax_nstruct=args.relax_nstruct,
                    seed=args.seed,
                )
                scan_df["source_pdb"] = os.path.basename(pdb_path)
                all_results.append(scan_df)
            except Exception as e:
                print(f"  ERROR on {pdb_path}: {e}")
                failed.append(sample)

        summary.append((antigen_chain, len(all_results), missing, failed))

        if not all_results:
            print(f"No results produced for chain {antigen_chain}.")
            continue

        combined = pd.concat(all_results, ignore_index=True)
        pivot = combined.pivot_table(
            index="source_pdb",
            columns="mutation",
            values="ddG_REU",
        )
        col_order = sorted(pivot.columns, key=lambda m: int("".join(filter(str.isdigit, m))))
        pivot = pivot[col_order]
        pivot.index.name = "source_pdb"
        pivot.columns.name = None

        print(f"\n--- ddG per structure (rows) x residue (columns), chain {antigen_chain} ---")
        print(pivot.to_string())

        out_path = os.path.join(
            analysis_dir, f"alanine_scan_{name}_{antigen_chain}{chunk_tag}.csv")
        pivot.to_csv(out_path)
        print(f"\nResults written to: {out_path}")

    # ── 4. Summary of this chunk ─────────────
    print(f"\n{'='*60}")
    print(f"Summary{f' of chunk {args.chunk_index}' if chunk_tag else ''} (seed {args.seed})")
    for antigen_chain, n_scanned, missing, failed in summary:
        print(f"  Chain {antigen_chain}: {n_scanned}/{len(list_best)} structures scanned")
        if missing:
            print(f"    no PDB file: {missing}")
        if failed:
            print(f"    failed:      {failed}")
