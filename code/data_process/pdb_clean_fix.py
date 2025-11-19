#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/08
# @Author  : Qiang Huang
# @File    : pdbfixer.py
import json
import os
import pymol
from pymol import cmd
pymol.finish_launching(['pymol', '-c'])
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def clean_structure(input_file, output_file, target_chains=None):
    """remove water, ions, and keep the certain chains"""
    cmd.load(input_file, "raw_structure")
    cmd.remove("resn HOH")
    cmd.remove("not polymer")

    # keep target chain
    if target_chains:
        cmd.remove(f"not chain {'+'.join(target_chains)}")

    cmd.save(output_file, "raw_structure")
    cmd.delete("raw_structure")

def fix_pdb(input_file, output_file):
    '''
    add missing resi and atoms of input pdb file
    '''
    fixer = PDBFixer(filename=input_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    with open(output_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

if __name__ =='__main__':
    with open('../data/chain_map_dict_test.json', 'r') as f:
        chain_map_dict = json.load(f)
    pdb_ids = chain_map_dict.keys()
    output_dir1 = '../data/cleaned_struc_test'
    output_dir2 = '../data/fixed_struc_test'
    os.makedirs(output_dir1, exist_ok=True)
    os.makedirs(output_dir2, exist_ok=True)
    for pdb_id in pdb_ids:
        input_raw = f'../data/tcr_struc_test/{pdb_id}.cif'
        output_clean = f'{output_dir1}/{pdb_id}.pdb'
        output_fix = f'{output_dir2}/{pdb_id}.pdb'
        target_chains = list(chain_map_dict[pdb_id].values())
        clean_structure(input_raw, output_clean, target_chains=target_chains)
        fix_pdb(output_clean, output_fix)

    # change the chain order A/B
    cmd.load('../data/fixed_struc_test/9o62.pdb', '9o62')
    cmd.alter("chain A", "chain='temp'")
    cmd.alter("chain B", "chain='A'")
    cmd.alter("chain temp", "chain='B'")
    cmd.sort()
    cmd.save('../data/fixed_struc_test/9o62.pdb', "9o62")
    cmd.delete('all')
