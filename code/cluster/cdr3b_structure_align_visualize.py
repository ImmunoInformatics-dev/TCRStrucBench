#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/25
# @Author  : Qiang Huang
# @File    :
import pandas as pd
from pymol import cmd, util
import os
import time
import logging
from tqdm import tqdm
'''
    align the CDR3beta structures
'''
def find_median(numbers):
    sorted_numbers = sorted(numbers)
    length = len(sorted_numbers)
    return sorted_numbers[length // 2]

if __name__ == '__main__':
    os.chdir('../../')
    models = ['AF3', 'tfold']
    for m in models:
        model_dir = f"../result/deepair/foldseek_{m}_clu"
        for clu in os.listdir(model_dir):
            input_dir = f"{model_dir}/{clu}"
            output_dir = f"{model_dir}/{clu}_aligned"
            os.makedirs(output_dir, exist_ok=True)

            logging.basicConfig(filename=f'{output_dir}/structure_alignment.log',
                                level=logging.INFO,
                                format='%(asctime)s - %(levelname)s - %(message)s')
            # file names
            structure_files = [f for f in os.listdir(input_dir)
                               if f.endswith((".pdb", ".cif"))]

            if not structure_files:
                logging.error("there is no PDB/CIF file")
                exit()

            logging.info(f"find {len(structure_files)} structure files")

            # load all struc
            structures = []
            pdb_size = []
            start_time = time.time()
            for i, filename in enumerate(tqdm(structure_files, desc="loading...")):
                try:
                    obj_name = os.path.splitext(filename)[0]
                    filepath = os.path.join(input_dir, filename)
                    stats_file = os.stat(filepath)
                    pdb_size.append(stats_file.st_size)
                    #
                    cmd.load(filepath, obj_name)
                    structures.append(obj_name)
                    logging.info(f"loaded: {filename} as {obj_name}")
                except Exception as e:
                    logging.error(f"{filename} loaded failed: {str(e)}")

            # ref structure
            median_id = pdb_size.index(find_median(pdb_size))
            reference = structures[median_id]
            logging.info(f"ref struc is: {reference}")

            # alignment
            alignment_results = []

            for i, structure in enumerate(tqdm(structures, desc="aligning")):
                if structure == reference:
                    continue

                try:
                    # align to ref.
                    rmsd = cmd.align(structure, reference)
                    alignment_results.append((structure, rmsd[0], rmsd[1]))
                    logging.info(f"align {structure} to {reference}: RMSD = {rmsd[0]:.3f}")
                except Exception as e:
                    logging.error(f"align {structure} to: {str(e)}")

            rmsd_values = [r[1] for r in alignment_results]
            structure_files.pop(median_id)
            rmsd_values_df = pd.DataFrame(rmsd_values, index=structure_files, columns=['RMSD'])
            rmsd_values_df.to_csv(f'{output_dir}/rmsd_values_df_250925.csv')
            # print log
            print("\n=== Abstract ===")
            print(f"Inference sturc: {reference}")
            print(f"aligned number: {len(alignment_results)}")
            print("RMSD distribution:")
            print(f"  Min: {min(rmsd_values):.3f} Å")
            print(f"  Max: {max(rmsd_values):.3f} Å")
            print(f"  Mean: {sum(rmsd_values) / len(rmsd_values):.3f} Å")

            ## visualize
            try:
                # create
                aligned_group = "aligned_structures"
                cmd.create(aligned_group, " or ".join(structures))

                # style
                cmd.show("cartoon", aligned_group)
                cmd.set("cartoon_transparency", 0.7, aligned_group)

                # others with gray
                for i, structure in enumerate(structures[1:]):
                    color_name = f"color_{i}"
                    cmd.set_color(color_name, [i / len(structures), 0.5, 1.0 - i / len(structures)])
                    util.color_carbon(color_name, structure)
                util.color_carbon("magenta", reference)

                # background
                cmd.bg_color("white")
                cmd.orient(aligned_group)
                cmd.zoom(aligned_group, 5)  # zoom 5 fold

                # save
                session_path = os.path.join(output_dir, "aligned_structures.pse")
                cmd.save(session_path)
                logging.info(f"PyMOL result was saved at: {session_path}")
                logging.info(f"Finished! total time: {time.time() - start_time:.2f} seconds")
                print(f"\nvisiualization was saved at: {session_path}")


            except Exception as e:
                logging.error(f"vision failure: {str(e)}")
                print(f"error in progress: {str(e)}")
            finally:
                # clear
                cmd.delete("all")