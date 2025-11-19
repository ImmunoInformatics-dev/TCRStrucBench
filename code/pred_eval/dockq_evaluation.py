#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TRCStruc
# @Purpose : evaluation DockQ
# @Time    : 2025/6/9
# @Author  : Qiang Huang
# @File    :
import pymol, json, os
import pandas as pd
import numpy as np
from pymol import cmd
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
pymol.finish_launching(['pymol', '-c'])

def select_res(full_seq, vdomain_seq):
    end5 = vdomain_seq[-5:-1]
    end = full_seq.find(end5) + 5 + 1
    return f'1-{end}'

def clean_hydro_trim(input_file, Vab_seq, output_file, pdb_id):
    cmd.load(input_file, pdb_id)
    cmd.remove("hydro")
    seqA = ''.join(cmd.get_fastastr('Chain A').splitlines()[1:])
    selA = select_res(seqA, Vab_seq[0])
    seqB = ''.join(cmd.get_fastastr('Chain B').splitlines()[1:])
    selB = select_res(seqB, Vab_seq[1])
    cmd.select('Vdomain', f'chain A and resi {selA} or (chain B and resi {selB})')
    cmd.save(output_file, 'Vdomain')
    cmd.delete('all')

def cal_dockq(pred_file, cry_file, chain_map = {"A": "A", "B": "B"}):
    model = load_PDB(pred_file)
    native = load_PDB(cry_file)
    dockq = run_on_all_native_interfaces(model, native, chain_map=chain_map)
    extract_keys = ['DockQ', 'iRMSD', 'LRMSD', 'fnat', 'F1']
    out_dq = [dockq[0]['AB'].get(x) for x in extract_keys]
    return out_dq

if __name__ == '__main__':
    os.chdir(r'D:\\USER\\Research\\TCRsructure\\code')
    tcr_anno = pd.read_csv("../result/annotation/tcr_anno_df.csv")
    tcr_names = tcr_anno['Name']
    pdb_ids2 = [x.split('_')[0] for x in tcr_names]
    pdb_ids = np.unique(pdb_ids2)
    # dict
    tcr_seq_dic = {}
    for kk, vv in zip(pdb_ids2, tcr_anno['Vdomain']):
        if kk in tcr_seq_dic.keys():
            tcr_seq_dic[kk].append(vv)
        else:
            tcr_seq_dic.update({kk: [vv]})

    #
    with open('../data/chain_map_dict_test.json', 'r') as f:
        chain_map_dict = json.load(f)
    input_dir = '../data/fixed_struc_test'
    output_dir = '../data/nohydro_trim_test'
    os.makedirs(output_dir, exist_ok=True)
    #
    dq_af2 = {}
    dq_af3 = {}
    dq_tcm2 = {}
    dq_tcb2 = {}
    dq_esm2 = {}
    dq_tfold = {}
    for pdb_id in pdb_ids:
        input_file = f'{input_dir}/{pdb_id}.pdb'
        output_clean = f'{output_dir}/{pdb_id}.pdb'
        Vab_seq = tcr_seq_dic[pdb_id]
        if not os.path.exists(output_clean):
            clean_hydro_trim(input_file, Vab_seq, output_clean, pdb_id)

        # AF2
        pred_af2 = f'../result/AF2/test_paired/{pdb_id}.pdb'
        dq_af2[pdb_id] = cal_dockq(pred_af2, output_clean, chain_map={"A": "A", "B": "B"})
        # AF3
        pred_af3 = f'../result/AF3/test_paired/{pdb_id}.cif'
        dq_af3[pdb_id] = cal_dockq(pred_af3, output_clean, chain_map={"A": "A", "B": "B"})
        # TCRmodel2
        pred_tcm2 = f'../result/tcrmodel2/test_paired/{pdb_id}.pdb'
        dq_tcm2[pdb_id] = cal_dockq(pred_tcm2, output_clean, chain_map={"A": "A", "B": "B"})
        # ESMFold2
        pred_esm2 = f'../result/ESMFold2/test_paired/{pdb_id}.pdb'
        dq_esm2[pdb_id] = cal_dockq(pred_esm2, output_clean, chain_map={"A": "A", "B": "B"})
        # tfold-tcr
        pred_tfold = f'../result/tfold/test_paired/{pdb_id}_TCR.pdb'
        dq_tfold[pdb_id] = cal_dockq(pred_tfold, output_clean, chain_map={"A": "A", "B": "B"})

        # aggregate
        col_names = ['PDB_ID', 'DockQ', 'iRMSD', 'LRMSD', 'fnat', 'F1']

        af2_dockq_df = pd.DataFrame(dq_af2).T.reset_index()
        af2_dockq_df.columns = col_names

        af3_dockq_df = pd.DataFrame(dq_af3).T.reset_index()
        af3_dockq_df.columns = col_names

        tcrmodel2_dockq_df = pd.DataFrame(dq_tcm2).T.reset_index()
        tcrmodel2_dockq_df.columns = col_names

        esmfold2_dockq_df = pd.DataFrame(dq_esm2).T.reset_index()
        esmfold2_dockq_df.columns = col_names

        tfold_dockq_df = pd.DataFrame(dq_tfold).T.reset_index()
        tfold_dockq_df.columns = col_names

        with pd.ExcelWriter(f'../result/dockq/dockq_summary_test_paired.xlsx', engine='openpyxl') as writer:
            af2_dockq_df.to_excel(writer, sheet_name='AF2')
            af3_dockq_df.to_excel(writer, sheet_name='AF3')
            esmfold2_dockq_df.to_excel(writer, sheet_name='ESMFoldv1')
            tcrmodel2_dockq_df.to_excel(writer, sheet_name='TCRmodel2')
            tfold_dockq_df.to_excel(writer, sheet_name='tfold-TCR')