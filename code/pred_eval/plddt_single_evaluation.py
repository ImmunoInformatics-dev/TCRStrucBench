#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TRCStruc
# @Purpose : evaluation with different index
# @Time    : 2025/6/23
# @Author  : Qiang Huang
# @File    : 
import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser
def region_dict_product(df, key_col, seq_col, FWs=True, pair_mode=True):
    for colname in key_col:
        df[f'{colname}_r'] = [get_region_range(x, y) for x, y in zip(df[colname], df[seq_col])]
    if FWs:
        df['FWs_r'] = [[a, b, c, d] for a, b, c, d in zip(df['FW1_r'], df['FW2_r'], df['FW3_r'], df['FW4_r'])]
    if pair_mode:
        df['PDB_ID'] = [x.split('_')[0] for x in df['Name']]
    return df
#
def get_single_chain_plddt(structure_file, chain_id, range_dic):
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

    chain_data = {}
    b_factors = []
    for chain in structure.get_chains():
        for residue in chain:
            res_id = residue.id
            if res_id[0] == ' ':
                residue_num = res_id[1]
                atom = next(residue.get_iterator(), None)
                if atom:
                    b_factors.append(atom.bfactor)
    chain_data[chain_id] = np.mean(b_factors)
    #
    region_plddt = {}
    for k, v in range_dic.items():
        if len(v) == 2:
            region_plddt[k] = np.mean(b_factors[v[0]:v[1]])
        else:  # fw1+fw2+fw3+fw4
            fws = []
            for z in v:
                fws.extend(b_factors[z[0]:z[1]])
            region_plddt[k] = np.mean(fws)
            
    return chain_data, region_plddt

def get_region_range(part_seq, whole_seq):
    start = whole_seq.find(part_seq[:5])  # first match
    end = start + len(part_seq)
    if end > len(whole_seq):
        end = len(whole_seq)
    return [start, end]
#
def stack_chain_out(plddt_df, index_col, pdb_col, chain_col, new_colname):
    '''
    :param plddt_df: the curated df of rmsd for all single chains
    :param index_col: the pdb_chain col
    :param pdb_col: name of pdb_id col
    :param chain_col: name of chain id col
    :return:
    '''
    plddt_df[pdb_col] = plddt_df[index_col].map(lambda x: x.split('_')[0])
    plddt_df[chain_col] = ['A', 'B'] * 2
    plddt_df.set_index([pdb_col, chain_col], inplace=True)
    plddt_df_w = plddt_df.unstack(level=1, sort=False)
    plddt_df_w = plddt_df_w.iloc[:, 2:]
    plddt_df_w = plddt_df_w.reset_index()
    plddt_df_w.columns = new_colname
    return plddt_df_w

if __name__ == '__main__':
    os.chdir(r'D:\\USER\\Research\\TCRsructure\\code')
    tcr_info = pd.read_csv('../result/annotation/tcr_anno_df.csv')
    tcr_info = region_dict_product(tcr_info,
                                   key_col=['FW1', 'FW2', 'FW3', 'FW4', 'CDR1', 'CDR2', 'CDR3'],
                                   seq_col='Vdomain', FWs=True, pair_mode=False)
    pdb_chains = tcr_info['Name']
    #
    af2_plddt = {}
    af3_plddt = {}
    tcrmodel2_plddt = {}
    tfold_plddt = {}
    esmfold2_plddt = {}
    for idx, pdb_chain in enumerate(pdb_chains):
        pdb_id = pdb_chain.split('_')[0]
        tcr_sele = tcr_info[tcr_info['Name'] == pdb_chain][tcr_info.columns[tcr_info.columns.str.endswith('_r')]]
        if idx/2 == 0:
            chain_id = 'A'
            tcr_sele.index = ['A']
        else:
            chain_id = 'B'
            tcr_sele.index = ['B']
        range_dict = tcr_sele.T.to_dict()[chain_id]
        # AF2
        file_af2 = f'../result/AF2/test_single/{pdb_chain}.pdb'  
        chain_data_af2, region_plddt_af2 = get_single_chain_plddt(file_af2, chain_id, range_dict)
        af2_plddt[pdb_chain] = list(chain_data_af2.values()) + list(pd.DataFrame(region_plddt_af2, index=[0]).T[0])
        # AF3
        file_af3 = f'../result/AF3/test_single/{pdb_chain}.cif'  
        chain_data_af3, region_plddt_af3 = get_single_chain_plddt(file_af3, chain_id, range_dict)
        af3_plddt[pdb_chain] = list(chain_data_af3.values()) + list(pd.DataFrame(region_plddt_af3, index=[0]).T[0])
        # TCRmodel2
        file_tcrm2 = f'../result/tcrmodel2/test_single/{pdb_chain}.pdb'  
        chain_data_tcrm2, region_plddt_tcrm2 = get_single_chain_plddt(file_tcrm2, chain_id, range_dict)
        tcrmodel2_plddt[pdb_chain] = list(chain_data_tcrm2.values()) + list(pd.DataFrame(region_plddt_tcrm2, index=[0]).T[0])
        # tfold-TCR
        file_tfold = f'../result/tfold/test_single/{pdb_chain}_TCR.pdb'  
        chain_data_tfold, region_plddt_tfold = get_single_chain_plddt(file_tfold, chain_id, range_dict)
        tfold_plddt[pdb_chain] = list(chain_data_tfold.values()) + list(pd.DataFrame(region_plddt_tfold, index=[0]).T[0])
        # ESMFold2
        file_esm2 = f'../result/ESMFold2/test_single/{pdb_chain}.pdb'  
        chain_data_esm2, region_plddt_esm2 = get_single_chain_plddt(file_esm2, chain_id, range_dict)
        esmfold2_plddt[pdb_chain] = list(chain_data_esm2.values()) + list(pd.DataFrame(region_plddt_esm2, index=[0]).T[0])


    #
    col_names = ['PDB_ID', 'Valpha', 'Vbeta', 'A-FW1', 'B-FW1', 'A-FW2', 'B-FW2', 'A-FW3', 'B-FW3',
                 'A-FW4', 'B-FW4', 'A-CDR1', 'B-CDR1', 'A-CDR2', 'B-CDR2', 'A-CDR3', 'B-CDR3', 'A-FWs', 'B-FWs']
    af2_plddt_df = pd.DataFrame(af2_plddt).T.reset_index()
    af2_plddt_df = stack_chain_out(af2_plddt_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain', new_colname=col_names)

    af3_plddt_df = pd.DataFrame(af3_plddt).T.reset_index()
    af3_plddt_df = stack_chain_out(af3_plddt_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain', new_colname=col_names)

    tcrmodel2_plddt_df = pd.DataFrame(tcrmodel2_plddt).T.reset_index()
    tcrmodel2_plddt_df = stack_chain_out(tcrmodel2_plddt_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain', new_colname=col_names)

    tfold_plddt_df = pd.DataFrame(tfold_plddt).T.reset_index()
    tfold_plddt_df = stack_chain_out(tfold_plddt_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain', new_colname=col_names)

    esmfold2_plddt_df = pd.DataFrame(esmfold2_plddt).T.reset_index()
    esmfold2_plddt_df = stack_chain_out(esmfold2_plddt_df, index_col='index', pdb_col='PDB_ID', chain_col='Chain', new_colname=col_names)

    with pd.ExcelWriter('../result/plddt/plddt_single_summary_test_single.xlsx', engine='openpyxl') as writer:
        af2_plddt_df.to_excel(writer, sheet_name='AF2')
        tcrmodel2_plddt_df.to_excel(writer, sheet_name='TCRmodel2')
        esmfold2_plddt_df.to_excel(writer, sheet_name='ESMFoldv1')
        tfold_plddt_df.to_excel(writer, sheet_name='tfold-TCR')
        af3_plddt_df.to_excel(writer, sheet_name='AF3')
