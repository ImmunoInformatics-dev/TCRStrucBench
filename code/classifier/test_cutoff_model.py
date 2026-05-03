#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/23
# @Author  : Qiang Huang
# @File    : cutoff_model_for_single.py

import os
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, accuracy_score

def create_multi_bars(datas, group, labels, tick_step=3, legend=True,
                      group_gap=0.3, bar_gap=0, ax=None, titles=None):

    x = np.arange(len(labels)) * tick_step
    group_num = len(datas)
    group_width = tick_step - group_gap
    bar_span = group_width / group_num
    bar_width = bar_span - bar_gap

    if ax is None:
        ax = plt.gca()

    for index, y in enumerate(datas):
        ax.bar(x + index*bar_span, y, bar_width, label=group[index])
    ax.set_ylabel('Scores')
    ax.set_title(titles)
    ticks = x + (group_width - bar_span) / 2
    ax.set_xticks(ticks, labels, rotation=45, ha='right')
    if legend:
        ax.legend(bbox_to_anchor=(1.05, 0.5), loc=6, borderaxespad=0)
    else:
        ax.legend().remove()
    return ax

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
    CURRENT_DIR = Path(__file__).resolve().parent
    PROJECT_ROOT = CURRENT_DIR.parent.parent
    fig_size = (16, 8)
    current_date = datetime.now()
    datestamp = f'{str(current_date.year)[2:]}{current_date.month:02d}{current_date.day:02d}'
    models = ['tfold-TCR', 'AF3']
    modes = {'single_cdr3': ['A-CDR3', 'B-CDR3'],
             'paired_cdr3': ['A-CDR3', 'B-CDR3'],
             'single': ['Valpha', 'Vbeta'],
             'paired': ['Vab']
             }

    # load dict
    model_testdict = {}
    label_dict = {}
    pred_dict = {}
    for model in models:
        with open(PROJECT_ROOT / f'result/train_model/{model}_cutoff_model_cv_3_dict_250922.json', 'r') as file:
            cutoff_dict = json.load(file)
        outdir = PROJECT_ROOT / f'result/test_model/{model}'
        os.makedirs(outdir, exist_ok=True)

        mode_ls = []
        for mode, fea_ls in modes.items():
            mod = mode.split('_')[0]
            for feat in fea_ls:
                mode_name = mod + '_' + feat
                mode_ls.append(mode_name)
        #
        sen_ls = []
        spe_ls = []
        acc_ls = []
        for idx, mode in enumerate(mode_ls):
            mode = mode_ls[idx]
            model_mtx = pd.read_excel(PROJECT_ROOT / f'result/test_model/make_{model}_modelling_mtx_250922.xlsx',
                                      sheet_name=mode, index_col=0)
            X = model_mtx.iloc[:, :-1].values
            y = model_mtx.iloc[:, -1].values
            label_dict[mode] = y

            model_mtx['pred'] = test_pred(X, mode, cutoff_dict)
            pred_df = model_mtx.loc[:, ['label', 'pred']]
            model_testdict[mode] = pred_df
            pred_dict[mode] = ['●' if x == 1 else '○' for x in model_mtx['pred']]
            #
            acc_ls.append(accuracy_score(y, model_mtx['pred']))
            # cm
            cm_df = pd.DataFrame(confusion_matrix(y, model_mtx['pred']), index=['0', '1'], columns=['0', '1'])
            sensitivity = cm_df.iloc[1, 1] / cm_df.sum(axis='columns')[1]
            sen_ls.append(sensitivity)
            specificity = cm_df.iloc[0, 0] / cm_df.sum(axis='columns')[0]
            spe_ls.append(specificity)

        hist_df = pd.DataFrame({'Sensitivity': sen_ls,
                                'Specificity': spe_ls,
                                'Accuracy': acc_ls},
                               index=mode_ls).T
        hist_df.to_excel(outdir / f'{model}_sen_spe_summary.xlsx')

        # histgram
        create_multi_bars(hist_df.values,
                          hist_df.index.values,
                          hist_df.columns,
                          bar_gap=0, group_gap=0.5)
        plt.savefig(outdir / f'{model}_sen_spe_acc.pdf',
                    dpi=300, bbox_inches='tight')
        plt.close()
