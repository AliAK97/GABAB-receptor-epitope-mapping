"""
Sort and rename AlphaFold3 output, and convert it for PyRosetta.

This script separates three AlphaFold3 output file types into separate 
folders, and replaces the seed and sample names with a simple index.
Files keep the name of the complex, and every indexed file remains 
traceable to the original seed and sample names.

Finally, each .cif is converted to a .pdb file, the format read by PyRosetta.
"""

import os
import shutil
import argparse
from Bio import PDB


def separate_files(pdb_id, output_files_path, sorted_cif_folder, sorted_json_folder, sorted_pae_folder):
    """Copies the .cif, summary .json and pae .json files of one complex into their own folders."""
    counter_cif = 0
    counter_json = 0
    counter_pae = 0
    base_path = f'{output_files_path}/{pdb_id}'
    dest_path_cif = sorted_cif_folder
    dest_path_json = sorted_json_folder
    dest_path_pae = sorted_pae_folder
    for first_level_name in os.listdir(base_path):
        first_level_path = os.path.join(base_path, first_level_name)
        if os.path.isdir(first_level_path):
            for file in os.listdir(first_level_path):
                src_file = os.path.join(first_level_path, file)
                if file.endswith('.cif'):
                    os.makedirs(f'{dest_path_cif}/{pdb_id}_cif', exist_ok=True)
                    dst_file = os.path.join(f'{dest_path_cif}/{pdb_id}_cif', file)
                    shutil.copy2(src_file, dst_file)
                    counter_cif += 1

                elif file.endswith('summary_confidences.json'):
                    os.makedirs(f'{dest_path_json}/{pdb_id}_json_summary', exist_ok=True)
                    dst_file = os.path.join(f'{dest_path_json}/{pdb_id}_json_summary', file)
                    shutil.copy2(src_file, dst_file)
                    counter_json += 1

                elif file.endswith('confidences.json') and 'summary' not in file:
                    os.makedirs(f'{dest_path_pae}/{pdb_id}_json_pae', exist_ok=True)
                    dst_file = os.path.join(f'{dest_path_pae}/{pdb_id}_json_pae', file)
                    shutil.copy2(src_file, dst_file)
                    counter_pae += 1

    print(f"Total .cif files copied: {counter_cif}")
    print(f"Total summary.json files copied: {counter_json}")
    print(f"Total pae.json files copied: {counter_pae}")


def rename_files(pdb_id, sorted_cif_folder, sorted_json_folder, sorted_pae_folder):
    """Renames the files in each folder to <pdb_id>_<index>, in sorted order."""
    cif_path  = f'{sorted_cif_folder}/{pdb_id}_cif'
    json_path = f'{sorted_json_folder}/{pdb_id}_json_summary'
    pae_path  = f'{sorted_pae_folder}/{pdb_id}_json_pae'

    try:
        for i, file in enumerate(sorted(os.listdir(cif_path))):
            if file.endswith('.cif'):
                new_name = f'{pdb_id}_{i}.cif'
                os.rename(os.path.join(cif_path, file), os.path.join(cif_path, new_name))

        for i, file in enumerate(sorted(os.listdir(json_path))):
            if file.endswith('.json'):
                new_name = f'{pdb_id}_{i}.json'
                os.rename(os.path.join(json_path, file), os.path.join(json_path, new_name))

        for i, file in enumerate(sorted(os.listdir(pae_path))):
            if file.endswith('.json'):
                new_name = f'{pdb_id}_{i}_pae.json'
                os.rename(os.path.join(pae_path, file), os.path.join(pae_path, new_name))
    except FileNotFoundError as e:
        print(f"WARNING: {e.filename} not found, renaming skipped for {pdb_id}")


if __name__ == "__main__":
    # ── 1. Arguments and paths ─────────────
    parser = argparse.ArgumentParser(
        description="Separate and rename AF3 output files."
    )
    parser.add_argument("structure_name",
                        help="Path to the structure directory")
    args = parser.parse_args()

    base = f"{str(args.structure_name).rstrip('/')}/full_antigen"
    output_files_path       = f"{base}/output"
    sorted_cif_folder       = f"{base}/sorted_output/cif_folder"
    sorted_json_folder      = f"{base}/sorted_output/json_folder"
    sorted_pae_folder       = f"{base}/sorted_output/pae_folder"
    sorted_pdb_folder       = f"{base}/sorted_output/pdb_folder"

    for folder in (sorted_cif_folder, sorted_json_folder, sorted_pae_folder, sorted_pdb_folder):
        os.makedirs(folder, exist_ok=True)

    # ── 2. Separate and rename each AF3 output folder ─────────────
    pae_files = os.listdir(sorted_pae_folder)
    pae_ids = {f.removesuffix('_json_pae').lower() for f in pae_files if not f.startswith('.')}

    for file in os.listdir(output_files_path):
        if file.startswith('.'):
            continue
        pdb_id = file
        if pdb_id.lower() in pae_ids:
            print(f"file {file} already exists")
        else:
            print(f"Processing PDB ID: {pdb_id}")
            separate_files(pdb_id, output_files_path, sorted_cif_folder, sorted_json_folder, sorted_pae_folder)
            rename_files(pdb_id, sorted_cif_folder, sorted_json_folder, sorted_pae_folder)
            print(f"==============Finished processing PDB ID: {pdb_id}============\n")

    # ── 3. Convert .cif to .pdb ─────────────
    mmcif_parser = PDB.MMCIFParser(QUIET=True, auth_residues=False)
    pdb_io = PDB.PDBIO()
    counter = 0

    for file_line_1 in os.listdir(sorted_cif_folder):
        first_line_path = os.path.join(sorted_cif_folder, file_line_1)

        if not os.path.isdir(first_line_path):
            continue

        output_subdir = os.path.join(sorted_pdb_folder, file_line_1.replace("_cif", "_pdb"))

        if os.path.exists(output_subdir):
            print(f"{file_line_1} exists")
            continue

        os.makedirs(output_subdir, exist_ok=True)

        for file_line_2 in os.listdir(first_line_path):
            second_line_path = os.path.join(first_line_path, file_line_2)

            if not file_line_2.endswith(".cif"):
                continue

            pdb_filename = file_line_2.replace(".cif", ".pdb")
            output_file_path = os.path.join(output_subdir, pdb_filename)

            structure = mmcif_parser.get_structure(file_line_2, second_line_path)
            pdb_io.set_structure(structure)
            pdb_io.save(output_file_path)

        print(f"{file_line_1} was turned into pdb")
        counter += 1

    print(f"Total converted structures: {counter}")
