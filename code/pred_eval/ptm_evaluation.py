#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    extract the ptm and iptm value at paired mode
'''
import os

import pandas as pd
import numpy as np
import json
from pathlib import Path
if __name__ == '__main__':
    CURRENT_DIR = Path(__file__).resolve().parent
    PROJECT_ROOT = CURRENT_DIR.parent.parent
    tcr_info = pd.read_csv(PROJECT_ROOT / 'result/annotation/tcr_anno_df.csv')
    pdb_ids = pd.unique(tcr_info['Name'].map(lambda x: x.split('_')[0]))
    ptm_dict = {'pdb_id': pdb_ids, 'AF2': [], 'AF3': [], 'TCRmodel2': [], 'ESMFoldv1': [], 'tfold-TCR': []}
    for pdb_id in pdb_ids:
        #
        af2_file = PROJECT_ROOT / f'result/AF2/test_paired/{pdb_id}_ranking.json'
        with open(af2_file, 'r') as f:
            af2_ptm_dict = json.load(f)
            ptm_dict['AF2'].append(af2_ptm_dict['iptm+ptm'][af2_ptm_dict['order'][0]])
        #
        af3_file = PROJECT_ROOT / f'result/AF3/test_paired/{pdb_id}_summary_confidences.json'
        with open(af3_file, 'r') as f:
            af3_ptm_dict = json.load(f)
            ptm_dict['AF3'].append(af3_ptm_dict['ranking_score'])
        #
        tcrm2_file = PROJECT_ROOT / f'result/tcrmodel2/test_paired/{pdb_id}_statistics.json'
        with open(tcrm2_file, 'r') as f:
            tcrm2_ptm_dict = json.load(f)
            if 'ranked_0' in tcrm2_ptm_dict.keys():
                ptm_dict['TCRmodel2'].append(tcrm2_ptm_dict['ranked_0']['model_confidence'])
            elif 'iptm+ptm' in tcrm2_ptm_dict.keys():
                ptm_dict['TCRmodel2'].append(tcrm2_ptm_dict['iptm+ptm'][tcrm2_ptm_dict['order'][0]])
            else:
                raise ValueError('check your input ptm files')
        #
        tfold_file = PROJECT_ROOT / f'result/tfold/test_paired/{pdb_id}_TCR.pdb'
        with open(tfold_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'Predicted pTM' in line:
                ptm_sc = float(line.split(' ')[-1])
            if 'Predicted ipTM' in line:
                iptm_sc = float(line.split(' ')[-1])
        if ptm_sc and iptm_sc:
            ptm_dict['tfold-TCR'].append(ptm_sc*0.2 + iptm_sc*0.8)
        else:
            ptm_dict['tfold-TCR'].append(None)
        #
        esm_ptm_dict = np.load(PROJECT_ROOT / f'result/ESMFold/test_paired/{pdb_id}.npy')
        ptm_dict['ESMFoldv1'].extend(list(esm_ptm_dict))
    ptm_df = pd.DataFrame(ptm_dict)

    with pd.ExcelWriter(PROJECT_ROOT / 'result/ptm/ptm_score_test_paired.xlsx', mode='w', engine='openpyxl') as writer:
        ptm_df.to_excel(writer, sheet_name='confidence score')
