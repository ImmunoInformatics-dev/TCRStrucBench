#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TRCStruc
# @Purpose : evaluation with different index
# @Time    : 2025/8/04
# @Author  : Qiang Huang
# @File    :
import os
import pymol
from pymol import cmd
import pandas as pd
import json

pymol.finish_launching(['pymol', '-c'])
def calculate_monomer_rmsd(pred_path, crystal_file, chain_id, domain_sel=None, backbone=False):
    """
    - chain_id: {'pred_chain': 'crystal_chain'}
    - domain_sel: {'chain': PyMOL GRAMMER, ...}
      EG.：{'A': 'resi 50-100', 'B': 'resi 30-80'}
    """
    pred_name = os.path.basename(pred_path).split('.')[0] + '_pred'
    crystal_name = os.path.basename(crystal_file).split('.')[0] + '_cry'

    cmd.load(pred_path, pred_name)
    cmd.load(crystal_file, crystal_name)

    if backbone:
        pred_sel = f"{pred_name} and name N+C+CA+O"
        crystal_sel = f"{crystal_name} and name N+C+CA+O"
    else:
        pred_sel = f"{pred_name}"
        crystal_sel = f"{crystal_name}"

    chain_rmsd = cmd.align(pred_sel, crystal_sel)[0]


    domain_rmsd = {}
    if domain_sel[chain_id]:
        if backbone:
            crystal_dom_sel = f"{crystal_name} and name N+C+CA+O"
        else:
            crystal_dom_sel = f"{crystal_name}"
        for region, ranges in domain_sel[chain_id].items():
            if type(ranges) == str:
                if backbone:
                    pred_dom_sel = f"{pred_name} and resi {ranges} and name N+C+CA+O"
                else:
                    pred_dom_sel = f"{pred_name} and resi {ranges}"
            else:  # fws
                if backbone:
                    pred_dom_sel = f"{pred_name} and (resi {ranges[0]} or resi {ranges[1]} or resi {ranges[2]} or resi {ranges[3]}) and name N+C+CA+O"
                else:
                    pred_dom_sel = f"{pred_name} and (resi {ranges[0]} or resi {ranges[1]} or resi {ranges[2]} or resi {ranges[3]})"

            if cmd.count_atoms(pred_dom_sel) == 0 or cmd.count_atoms(crystal_dom_sel) == 0:
                domain_rmsd[f'Region {region}'] = ''
                continue

            domain_rmsd[f'Region {region}'] = cmd.align(pred_dom_sel, crystal_dom_sel)[0]

    cmd.delete(pred_name)
    cmd.delete(crystal_name)

    return chain_rmsd, domain_rmsd

def get_region_range(part_seq, whole_seq):
    '''
    return the resi number of domains
    '''
    start = whole_seq.find(part_seq[:5])  # first match
    end = start + len(part_seq)
    start += 1
    if end > len(whole_seq):
        end = len(whole_seq)
    return f'{start}-{end}'

def stack_chain_out(rmsd_df, index_col, pdb_col, chain_col, new_colname):
    '''
    :param rmsd_df: the curated df of rmsd for all single chains
    :param index_col: the pdb_chain col
    :param pdb_col: name of pdb_id col
    :param chain_col: name of chain id col
    :return:
    '''
    rmsd_df[pdb_col] = rmsd_df[index_col].map(lambda x: x.split('_')[0])
    rmsd_df[chain_col] = ['A', 'B'] * 2
    rmsd_df.set_index([pdb_col, chain_col], inplace=True)
    rmsd_df_w = rmsd_df.unstack(level=1, sort=False)
    rmsd_df_w = rmsd_df_w.iloc[:, 2:]
    rmsd_df_w = rmsd_df_w.reset_index()
    rmsd_df_w.columns = new_colname
    return rmsd_df_w

if __name__ == '__main__':
    #
    os.chdir(r'D:\\USER\\Research\\TCRsructure\\code')
    tcr_info = pd.read_csv('../result/annotation/tcr_anno_df.csv')
    tcr_info['FW1_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW1'], tcr_info['Vdomain'])]
    tcr_info['FW2_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW2'], tcr_info['Vdomain'])]
    tcr_info['FW3_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW3'], tcr_info['Vdomain'])]
    tcr_info['FW4_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW4'], tcr_info['Vdomain'])]
    tcr_info['FWs_r'] = [[a, b, c, d] for a, b, c, d in zip(tcr_info['FW1_r'], tcr_info['FW2_r'], tcr_info['FW3_r'], tcr_info['FW4_r'])]
    tcr_info['CDR1_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['CDR1'], tcr_info['Vdomain'])]
    tcr_info['CDR2_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['CDR2'], tcr_info['Vdomain'])]
    tcr_info['CDR3_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['CDR3'], tcr_info['Vdomain'])]
    pdb_chains = tcr_info['Name']

    fixed_dir = '../data/fixed_struc'

    af2_rmsd = {}
    af3_rmsd = {}
    tcrmodel2_rmsd = {}
    tfold_rmsd = {}
    esmfold2_rmsd = {}
    for idx, pdb_chain in enumerate(pdb_chains):
        pdb_id = pdb_chain.split('_')[0]
        crystal_file = f'{fixed_dir}/{pdb_id}.pdb'
        tcr_sele = tcr_info[tcr_info['Name'] == pdb_chain][tcr_info.columns[tcr_info.columns.str.endswith('_r')]]
        if idx/2 == 0:
            chain_id = 'A'
            tcr_sele.index = ['A']
        else:
            chain_id = 'B'
            tcr_sele.index = ['B']
        range_dict = tcr_sele.T.to_dict()
        # af2
        af2_file = f'../result/AF2/test_single/{pdb_chain}.pdb'
        chain_rmsd_af2, domain_rmsd_af2 = calculate_monomer_rmsd(af2_file, crystal_file, chain_id, range_dict)
        af2_rmsd[pdb_chain] = [chain_rmsd_af2] + list(domain_rmsd_af2.values())
        # af3
        af3_file = f'../result/AF3/test_single/{pdb_chain}.cif'
        chain_rmsd_af3, domain_rmsd_af3 = calculate_monomer_rmsd(af3_file, crystal_file, chain_id, range_dict)
        af3_rmsd[pdb_chain] = [chain_rmsd_af3] + list(domain_rmsd_af3.values())
        # TCRmodel2
        tcrm2_file = f'../result/tcrmodel2/test_single/{pdb_chain}.pdb'
        chain_rmsd_tcrm2, domain_rmsd_tcrm2 = calculate_monomer_rmsd(tcrm2_file, crystal_file, chain_id, range_dict)
        tcrmodel2_rmsd[pdb_chain] = [chain_rmsd_tcrm2] + list(domain_rmsd_tcrm2.values())

        # ESMFold2
        esmfold2_file = f'../result/ESMFold2/test_single/{pdb_chain}.pdb'
        chain_rmsd_esm2, domain_rmsd_esm2 = calculate_monomer_rmsd(esmfold2_file, crystal_file, chain_id, range_dict)
        esmfold2_rmsd[pdb_chain] = [chain_rmsd_esm2] + list(domain_rmsd_esm2.values())

        # tfold-TCR
        tfold_file = f'../result/tfold/test_single/{pdb_chain}_TCR.pdb'
        chain_rmsd_tfold, domain_rmsd_tfold = calculate_monomer_rmsd(tfold_file, crystal_file, chain_id, range_dict)
        tfold_rmsd[pdb_chain] = [chain_rmsd_tfold] + list(domain_rmsd_tfold.values())

    #
    col_names = ['PDB_ID', 'Valpha', 'Vbeta', 'A-FW1', 'B-FW1', 'A-FW2', 'B-FW2', 'A-FW3', 'B-FW3',
                 'A-FW4', 'B-FW4', 'A-FWs', 'B-FWs', 'A-CDR1', 'B-CDR1', 'A-CDR2', 'B-CDR2', 'A-CDR3', 'B-CDR3']
    af2_rmsd_df = pd.DataFrame(af2_rmsd).T.reset_index()
    af2_rmsd_df = stack_chain_out(af2_rmsd_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain', new_colname=col_names)

    af3_rmsd_df = pd.DataFrame(af3_rmsd).T.reset_index()
    af3_rmsd_df = stack_chain_out(af3_rmsd_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain',
                                  new_colname=col_names)

    tcrmodel2_rmsd_df = pd.DataFrame(tcrmodel2_rmsd).T.reset_index()
    tcrmodel2_rmsd_df = stack_chain_out(tcrmodel2_rmsd_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain', new_colname=col_names)

    tfold_rmsd_df = pd.DataFrame(tfold_rmsd).T.reset_index()
    tfold_rmsd_df = stack_chain_out(tfold_rmsd_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain',
                                  new_colname=col_names)

    esmfold2_rmsd_df = pd.DataFrame(esmfold2_rmsd).T.reset_index()
    esmfold2_rmsd_df = stack_chain_out(esmfold2_rmsd_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain',
                                    new_colname=col_names)

    with pd.ExcelWriter(f'../result/rmsd/rmsd_summary_test_single_0804.xlsx', engine='openpyxl') as writer:
        af2_rmsd_df.to_excel(writer, sheet_name='AF2')
        af3_rmsd_df.to_excel(writer, sheet_name='AF3')
        esmfold2_rmsd_df.to_excel(writer, sheet_name='ESMFoldv1')
        tcrmodel2_rmsd_df.to_excel(writer, sheet_name='TCRmodel2')
        tfold_rmsd_df.to_excel(writer, sheet_name='tfold-TCR')