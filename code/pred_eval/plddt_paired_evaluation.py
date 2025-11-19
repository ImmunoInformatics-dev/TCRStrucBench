#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TRCStruc
# @Purpose : evaluation with different index
# @Time    : 2025/8/04
# @Author  : Qiang Huang
# @File    : 
import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser

#
def get_pair_chain_plddt(structure_file, range_dic):
    """
    :param structure_file: the structure file
    :param range_dic: the CDRs and FWs range for two chains, eg. {'A': {'Vdomain': [0, 126], 'CDR1': [27, 39]}, 'B': {},}
    :return: plddt_data (plddts dict of all chains and all residues),
    """
    if structure_file.endswith('.cif'):
        parser = MMCIFParser()
    else:
        parser = PDBParser()
    structure = parser.get_structure('model', structure_file)

    plddt_data = {}
    chain_data = {}
    for chain in structure.get_chains():
        chain_id = chain.id
        b_factors = []
        for residue in chain:
            res_id = residue.id
            if res_id[0] == ' ':
                residue_num = res_id[1]
                atom = next(residue.get_iterator(), None)
                if atom:
                    b_factors.append(atom.bfactor)
        plddt_data[chain_id] = b_factors
        chain_data[chain_id] = np.mean(b_factors)
    #
    avg_ab = np.mean(list(plddt_data.values())[0] + list(plddt_data.values())[1])

    region_plddt = {'A': {}, 'B': {}}
    for k, v in range_dic['A'].items():
        if len(v) == 2:
            region_plddt['A'][k] = np.mean(plddt_data['A'][v[0]:v[1]])
        else:  # fw1+fw2+fw3+fw4
            fws = []
            for z in v:
                fws.extend(plddt_data['A'][z[0]:z[1]])
            region_plddt['A'][k] = np.mean(fws)

    for k, v in range_dic['B'].items():
        if len(v) == 2:
            region_plddt['B'][k] = np.mean(plddt_data['B'][v[0]:v[1]])
        else:  # fw1+fw2+fw3+fw4
            fws = []
            for z in v:
                fws.extend(plddt_data['B'][z[0]:z[1]])
            region_plddt['B'][k] = np.mean(fws)
            
    return plddt_data, avg_ab, chain_data, region_plddt

def get_region_range(part_seq, whole_seq):
    start = whole_seq.find(part_seq[:5])  # first match
    end = start + len(part_seq)
    if end > len(whole_seq):
        end = len(whole_seq)
    return [start, end]
#
if __name__ == '__main__':
    os.chdir(r'D:\\USER\\Research\\TCRsructure\\code')
    tcr_info = pd.read_csv('../result/annotation/tcr_anno_df.csv')
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
    #
    af2_plddt = {}
    af3_plddt = {}
    tcrmodel2_plddt = {}
    tfold_plddt = {}
    esmfold2_plddt = {}
    for pdb_id in pdb_ids:
        tcr_sele = tcr_info[tcr_info['PDB_ID'] == pdb_id][tcr_info.columns[tcr_info.columns.str.endswith('_r')]]
        tcr_sele.index = ['A', 'B']
        range_dict = tcr_sele.T.to_dict()
        # AF2
        file_af2 = f'../result/AF2/test_paired/{pdb_id}.pdb' 
        _, avg_ab_af2, chain_data_af2, region_plddt_af2 = get_pair_chain_plddt(file_af2, range_dict)
        af2_plddt[pdb_id] = [avg_ab_af2] + list(chain_data_af2.values()) + list(pd.DataFrame(region_plddt_af2).unstack(level=1).values)
        # AF3
        file_af3 = f'../result/AF3/test_paired/{pdb_id}.cif' 
        _, avg_ab_af3, chain_data_af3, region_plddt_af3 = get_pair_chain_plddt(file_af3, range_dict)
        af3_plddt[pdb_id] = [avg_ab_af3] + list(chain_data_af3.values()) + list(pd.DataFrame(region_plddt_af3).unstack(level=1).values)
        # TCRmodel2
        file_tcrm2 = f'../result/tcrmodel2/test_paired/{pdb_id}.pdb' 
        _, avg_ab_tcrm2, chain_data_tcrm2, region_plddt_tcrm2 = get_pair_chain_plddt(file_tcrm2, range_dict)
        tcrmodel2_plddt[pdb_id] = [avg_ab_tcrm2] + list(chain_data_tcrm2.values()) + list(pd.DataFrame(region_plddt_tcrm2).unstack(level=1).values)
        # tfold-TCR
        file_tfold = f'../result/tfold/test_paired/{pdb_id}_TCR.pdb' 
        _, avg_ab_tfold, chain_data_tfold, region_plddt_tfold = get_pair_chain_plddt(file_tfold, range_dict)
        tfold_plddt[pdb_id] = [avg_ab_tfold] + list(chain_data_tfold.values()) + list(pd.DataFrame(region_plddt_tfold).unstack(level=1).values)
        # ESMFold2
        file_esm2 = f'../result/ESMFold2/test_paired/{pdb_id}.pdb' 
        _, avg_ab_esm2, chain_data_esm2, region_plddt_esm2 = get_pair_chain_plddt(file_esm2, range_dict)
        esmfold2_plddt[pdb_id] = [avg_ab_esm2] + list(chain_data_esm2.values()) + list(pd.DataFrame(region_plddt_esm2).unstack(level=1).values)


    #
    col_names = ['PDB_ID', 'Vab', 'Valpha', 'Vbeta', 'A-FW1', 'A-FW2', 'A-FW3', 'A-FW4', 'A-FWs', 'A-CDR1', 'A-CDR2', 'A-CDR3',
                 'B-FW1', 'B-FW2', 'B-FW3', 'B-FW4', 'B-FWs', 'B-CDR1', 'B-CDR2', 'B-CDR3']
    af2_plddt_df = pd.DataFrame(af2_plddt).T.reset_index()
    af2_plddt_df.columns = col_names

    af3_plddt_df = pd.DataFrame(af3_plddt).T.reset_index()
    af3_plddt_df.columns = col_names

    tcrmodel2_plddt_df = pd.DataFrame(tcrmodel2_plddt).T.reset_index()
    tcrmodel2_plddt_df.columns = col_names

    tfold_plddt_df = pd.DataFrame(tfold_plddt).T.reset_index()
    tfold_plddt_df.columns = col_names

    esmfold2_plddt_df = pd.DataFrame(esmfold2_plddt).T.reset_index()
    esmfold2_plddt_df.columns = col_names

    with pd.ExcelWriter('../result/plddt/plddt_summary_test_paired.xlsx', engine='openpyxl') as writer:
        af2_plddt_df.to_excel(writer, sheet_name='AF2')
        tcrmodel2_plddt_df.to_excel(writer, sheet_name='TCRmodel2')
        esmfold2_plddt_df.to_excel(writer, sheet_name='ESMFoldv1')
        tfold_plddt_df.to_excel(writer, sheet_name='tfold-TCR')
        af3_plddt_df.to_excel(writer, sheet_name='AF3')
