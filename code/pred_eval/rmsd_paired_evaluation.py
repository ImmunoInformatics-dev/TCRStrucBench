#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TRCStruc
# @Purpose : evaluation with different index
# @Time    : 2025/9/8
# @Author  : Qiang Huang
# @File    : rmsd_evaluation_250422.py
import os
import pymol
from pymol import cmd
import pandas as pd
import json

pymol.finish_launching(['pymol', '-c'])

def calculate_dimer_rmsd(pred_path, crystal_file, chain_map, domain_sel=None, backbone=False):
    """
    - chain_map: {'pred_chain': 'crystal_chain'}
    - domain_sel: {'chain': PyMOL选择语法, ...}
      示例：{'A': 'resi 50-100', 'B': 'resi 30-80'}
    """
    results = []
    pred_name = os.path.basename(pred_path).split('.')[0] + '_pred'
    crystal_name = os.path.basename(crystal_file).split('.')[0] + '_cry'

    # 加载结构
    cmd.load(pred_path, pred_name)
    cmd.load(crystal_file, crystal_name)

    # 分链对齐基础结构
    chain_rmsd = {}
    for p_chain, c_chain in chain_map.items():
        if backbone:
            pred_sel = f"{pred_name} and chain {p_chain} and name N+C+CA+O"
            crystal_sel = f"{crystal_name} and chain {c_chain} and name N+C+CA+O"
        else:
            pred_sel = f"{pred_name} and chain {p_chain}"
            crystal_sel = f"{crystal_name} and chain {c_chain}"
        chain_rmsd[f'Chain {p_chain}'] = cmd.align(pred_sel, crystal_sel)[0]

    # 整体复合体RMSD
    full_rmsd = cmd.align(pred_name, crystal_name)[0]

    # 结构域RMSD（需独立对齐）
    domain_rmsd = {'Chain A': {}, 'Chain B': {}}
    if domain_sel:
        for chain in domain_sel.keys():  # 分链
            if backbone:
                crystal_dom_sel = f"{crystal_name} and chain {chain_map[chain]} and name N+C+CA+O"
            else:
                crystal_dom_sel = f"{crystal_name} and chain {chain_map[chain]}"
            for region, ranges in domain_sel[chain].items():  # 分区域
                # 生成选择表达式
                if type(ranges) == str:
                    if backbone:
                        pred_dom_sel = f"{pred_name} and chain {chain} and resi {ranges} and name N+C+CA+O"
                    else:
                        pred_dom_sel = f"{pred_name} and chain {chain} and resi {ranges}"
                else:  # fws
                    if backbone:
                        pred_dom_sel = f"{pred_name} and chain {chain} and (resi {ranges[0]} or resi {ranges[1]} or resi {ranges[2]} or resi {ranges[3]}) and name N+C+CA+O"
                    else:
                        pred_dom_sel = f"{pred_name} and chain {chain} and (resi {ranges[0]} or resi {ranges[1]} or resi {ranges[2]} or resi {ranges[3]})"
                # crystal_dom_sel = f"{crystal_name} and chain {chain_map[chain]} and {sel_syntax}"


                # 验证选择是否有效
                if cmd.count_atoms(pred_dom_sel) == 0 or cmd.count_atoms(crystal_dom_sel) == 0:
                    domain_rmsd[f'Chain {chain}'][f'Region {region}'] = ''
                    continue

                # 独立对齐结构域
                domain_rmsd[f'Chain {chain}'][f'Region {region}'] = cmd.align(pred_dom_sel, crystal_dom_sel)[0]

    # 记录结果
    results.append({
        'file': pred_name,
        'full_rmsd': full_rmsd,
        'chain_rmsd': chain_rmsd,
        'domain_rmsd': domain_rmsd
    })
    cmd.delete(pred_name)
    cmd.delete(crystal_name)

    return full_rmsd, chain_rmsd, domain_rmsd

def get_region_range(part_seq, whole_seq):
    '''
    返回pdb结构中残基的位置
    '''
    start = whole_seq.find(part_seq[:5])  # first match
    end = start + len(part_seq)
    start += 1  # pdb中不使用切片
    if end > len(whole_seq):
        end = len(whole_seq)
    return f'{start}-{end}'

if __name__ == '__main__':
    #
    os.chdir(r'D:\\Z440\\Postdoctoral\\Research\\TCRsructure\\code')
    tcr_info = pd.read_csv('../result/annotation/tcr_anno_df_v2_250908.csv')
    tcr_info['PDB_ID'] = [x.split('_')[0] for x in tcr_info['Name']]
    tcr_info['FW1_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW1'], tcr_info['Vdomain'])]
    tcr_info['FW2_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW2'], tcr_info['Vdomain'])]
    tcr_info['FW3_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW3'], tcr_info['Vdomain'])]
    tcr_info['FW4_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['FW4'], tcr_info['Vdomain'])]
    tcr_info['FWs_r'] = [[a, b, c, d] for a, b, c, d in zip(tcr_info['FW1_r'], tcr_info['FW2_r'], tcr_info['FW3_r'], tcr_info['FW4_r'])]
    tcr_info['CDR1_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['CDR1'], tcr_info['Vdomain'])]
    tcr_info['CDR2_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['CDR2'], tcr_info['Vdomain'])]
    tcr_info['CDR3_r'] = [get_region_range(x, y) for x, y in zip(tcr_info['CDR3'], tcr_info['Vdomain'])]
    pdb_ids = tcr_info['PDB_ID'].unique()

    fixed_dir = '../data/fixed_struc_4'

    chain_map = {'A': 'A', 'B': 'B'}  # 经过pdb_clean_fixer后，晶体结构全变成双链AB了
    af2_rmsd = {}
    af3_rmsd = {}
    tcrmodel2_rmsd = {}
    tfold_rmsd = {}
    esmfold2_rmsd = {}
    for pdb_id in pdb_ids:
        # chain_map = chain_map_dict[pdb_id]
        crystal_file = f'{fixed_dir}/{pdb_id}.pdb'
        tcr_sele = tcr_info[tcr_info['PDB_ID'] == pdb_id][tcr_info.columns[tcr_info.columns.str.endswith('_r')]]
        tcr_sele.index = ['A', 'B']
        range_dict = tcr_sele.T.to_dict()
        # af3
        af3_file = f'../result/AF3/test2_paired/{pdb_id}.cif'
        full_rmsd_af3, chain_rmsd_af3, domain_rmsd_af3 = calculate_dimer_rmsd(af3_file, crystal_file, chain_map, range_dict)
        af3_rmsd[pdb_id] = [full_rmsd_af3] + list(chain_rmsd_af3.values()) + list(pd.DataFrame(domain_rmsd_af3).unstack(level=1).values)
        
        # tfold-TCR
        tfold_file = f'../result/tfold/test2_paired/{pdb_id}_TCR.pdb'
        full_rmsd_tfold, chain_rmsd_tfold, domain_rmsd_tfold = calculate_dimer_rmsd(tfold_file, crystal_file, chain_map, range_dict)
        tfold_rmsd[pdb_id] = [full_rmsd_tfold] + list(chain_rmsd_tfold.values()) + list(pd.DataFrame(domain_rmsd_tfold).unstack(level=1).values)

        #
        col_names = ['PDB_ID', 'Vab', 'Valpha', 'Vbeta', 'A-FW1', 'A-FW2', 'A-FW3', 'A-FW4', 'A-FWs', 'A-CDR1', 'A-CDR2', 'A-CDR3',
                     'B-FW1', 'B-FW2', 'B-FW3', 'B-FW4', 'B-FWs', 'B-CDR1', 'B-CDR2', 'B-CDR3']

        af2_rmsd_df = pd.DataFrame(af2_rmsd).T.reset_index()
        af2_rmsd_df.columns = col_names

        af3_rmsd_df = pd.DataFrame(af3_rmsd).T.reset_index()
        af3_rmsd_df.columns = col_names

        tcrmodel2_rmsd_df = pd.DataFrame(tcrmodel2_rmsd).T.reset_index()
        tcrmodel2_rmsd_df.columns = col_names

        tfold_rmsd_df = pd.DataFrame(tfold_rmsd).T.reset_index()
        tfold_rmsd_df.columns = col_names

        esmfold2_rmsd_df = pd.DataFrame(esmfold2_rmsd).T.reset_index()
        esmfold2_rmsd_df.columns = col_names

        with pd.ExcelWriter(f'../result/rmsd/rmsd_summary_test2_paired_0908.xlsx', engine='openpyxl') as writer:
            af3_rmsd_df.to_excel(writer, sheet_name='AF3')
            tfold_rmsd_df.to_excel(writer, sheet_name='tfold-TCR')