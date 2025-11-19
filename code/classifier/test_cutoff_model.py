#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/25
# @Author  : Qiang Huang
# @File    :
'''
'''
import os
from datetime import datetime
import numpy as np
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
def test_pred(X_data, mode, cutoff_dict):
    mod = mode.split('_')[0]
    predictions = np.zeros(len(X_data))
    cutoffs = cutoff_dict[mode]["best_cutoff_combination"]
    if 'CDR3' in mode:
        condition = (X_data[:, 0] > cutoffs[0]) & (X_data[:, 1] < cutoffs[1])
    elif mod == 'paired':
        condition = (X_data[:, 0] > cutoffs[0]) & (X_data[:, 1] > cutoffs[1])
    else:
        condition = X_data[:, 0] > cutoffs
    predictions[condition] = 1
    return predictions

if __name__ == "__main__":
    os.chdir('../../')
    fig_size = (16, 8)
    current_date = datetime.now()
    datestamp = f'{current_date.year}{current_date.month:02d}{current_date.day:02d}'
    models = ['tfold', 'AF3']
    mode = 'single_B-CDR3'

    # load dict
    model_testdict = {}
    label_dict = {}
    pred_dict = {}
    for model in models:
        if model == 'tfold':
            model2 = 'tfold-TCR'
        else:
            model2 = model
        with open(f'../result/train_model/simplify/{model2}_cutoff_model_cv_3_dict.json', 'r') as file:
            cutoff_dict = json.load(file)
        outdir = f'../result/test_model/simplify/deepair_{model}'
        os.makedirs(outdir, exist_ok=True)

        model_mtx = pd.read_excel(f"../result/deepair/plddt/{model}_single_B_plddt.xlsx", index_col=0)
        X = model_mtx.values

        model_mtx['pred'] = test_pred(X, mode, cutoff_dict)
        #
        model_mtx.to_excel(f'{outdir}/deepair_{model}_cutoff_discrimination_{datestamp}.xlsx')


