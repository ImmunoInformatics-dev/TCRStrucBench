#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/22
# @Author  : Qiang Huang
# @File    :
import os
from datetime import datetime
import numpy as np
import pandas as pd
import json
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, roc_curve, auc
import itertools

np.random.seed(42)

def get_candidate_cutoffs(feature_series, n_quantiles=20):
    unique_vals = np.unique(feature_series)
    if len(unique_vals) > n_quantiles:
        candidate_cutoffs = np.unique(np.quantile(feature_series, np.linspace(0, 1, n_quantiles)))
    else:
        candidate_cutoffs = (unique_vals[:-1] + unique_vals[1:]) / 2.0
    return candidate_cutoffs

def predict_with_cutoff(X_data, cutoff_f1):
    predictions = np.zeros(len(X_data))
    condition = (X_data[:, 0] > cutoff_f1)
    predictions[condition] = 1
    return predictions

def predict_with_2cutoffs(X_data, cutoff_f1, cutoff_f2, cdr3=False):
    predictions = np.zeros(len(X_data))
    if cdr3:
        condition = (X_data[:, 0] > cutoff_f1) & (X_data[:, 1] < cutoff_f2)
    else:
        condition = (X_data[:, 0] > cutoff_f1) & (X_data[:, 1] > cutoff_f2)
    predictions[condition] = 1
    return predictions

def predict_new(features, cutoffs, cdr3=False):
    if cdr3:
        label = int((features[0] > cutoffs[0]) & (features[1] < cutoffs[1]))
    elif len(features) > 1:
        label = int(all(x > y for x, y in zip(features, cutoffs)))
    else:
        label = int(features > cutoffs)
    return label
#
if __name__ == "__main__":
    os.chdir('../../')
    current_date = datetime.now()
    datestamp = f'{current_date.year}{current_date.month:02d}{current_date.day:02d}'
    models = ['tfold-TCR', 'AF3']
    modes = {'single_cdr3': ['A-CDR3', 'B-CDR3'],
             'paired_cdr3': ['A-CDR3', 'B-CDR3'],
             'paired': ['Vab'],
             'single': ['Valpha', 'Vbeta']
             }

    # save dict
    model_dict = {}
    for model in models:
        for mode, fea_ls in modes.items():
            mod = mode.split('_')[0]
            for feat in fea_ls:
                mode_name = mod + '_' + feat
                outdir = f'../result/cutoff_grid/{model}/{datestamp}'
                os.makedirs(outdir, exist_ok=True)
                model_mtx = pd.read_excel(f'../result/train_model/simplify/make_{model}_modelling_mtx_250922.xlsx',
                                          sheet_name=mode_name, index_col=0)
                X = model_mtx.iloc[:, :-1].values
                y = model_mtx.iloc[:, -1].values
                if mode == 'single':
                    cutoffs_f1 = get_candidate_cutoffs(X[:, 0])
                # elif 'cdr3' in mode:
                else:
                    cutoffs_f1 = get_candidate_cutoffs(X[:, 0])
                    cutoffs_f2 = get_candidate_cutoffs(X[:, 1])
                # data split
                n_splits = 3
                skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
                best_mean_auc = 0
                best_cutoff_combination = None
                cv3_results = None
                #
                if mode != 'single':
                    all_cutoff_combinations = list(itertools.product(cutoffs_f1, cutoffs_f2))
                    # 遍历cutoff组合
                    for cutoff_f1, cutoff_f2 in all_cutoff_combinations:
                        cv_scores = []
                        for train_index, val_index in skf.split(X, y):
                            X_train, X_val = X[train_index, :], X[val_index, :]
                            y_train, y_val = y[train_index], y[val_index]
                            # 区分
                            if 'cdr3' in mode:
                                y_val_pred = predict_with_2cutoffs(X_val, cutoff_f1, cutoff_f2, cdr3=True)
                            else:
                                y_val_pred = predict_with_2cutoffs(X_val, cutoff_f1, cutoff_f2)

                            # ROC_AUC
                            FPR, TPR, _ = roc_curve(y_val, y_val_pred)
                            score = auc(FPR, TPR)
                            cv_scores.append(score)

                        mean_score = np.mean(cv_scores)
                        # 更新
                        if mean_score > best_mean_auc:
                            best_mean_auc = mean_score
                            best_cutoff_combination = (cutoff_f1, cutoff_f2)
                            cv3_results = cv_scores
                else:
                    # cutoff
                    for cutoff_f1 in cutoffs_f1:
                        cv_scores = []
                        for train_index, val_index in skf.split(X, y):
                            X_train, X_val = X[train_index, :], X[val_index, :]
                            y_train, y_val = y[train_index], y[val_index]
                            #
                            y_val_pred = predict_with_cutoff(X_val, cutoff_f1)

                            # ROCAUC
                            FPR, TPR, _ = roc_curve(y_val, y_val_pred)
                            score = auc(FPR, TPR)
                            cv_scores.append(score)

                        mean_score = np.mean(cv_scores)  # 5折内平均值
                        #
                        if mean_score > best_mean_auc:
                            best_mean_auc = mean_score
                            best_cutoff_combination = cutoff_f1
                            cv3_results = cv_scores

                model_dict[mode_name] = {'best_mean_auc': best_mean_auc,
                                         'best_cutoff_combination': best_cutoff_combination,
                                         'cv3_results': cv3_results
                                         }
        with open(f'../result/train_model/simplify/{model}_cutoff_model_cv_{n_splits}_dict.json', mode='w', encoding='utf8') as f:
            json.dump(model_dict, f, ensure_ascii=False, indent=4)
        # excel
        best_mean_auc_ls = []
        cutoff_1_ls = []
        cutoff_2_ls = []

        for sh_name, eval_dict in model_dict.items():
            best_mean_auc_ls.append(eval_dict['best_mean_auc'])
            if 'single_V' in sh_name:
                cutoff_1_ls.append(eval_dict['best_cutoff_combination'])
                cutoff_2_ls.append(None)
            else:
                cutoff_1_ls.append(eval_dict['best_cutoff_combination'][0])
                cutoff_2_ls.append(eval_dict['best_cutoff_combination'][1])
        summary_df = pd.DataFrame({'Modes': model_dict.keys(),
                                   'Mean_AUC': best_mean_auc_ls,
                                   'Cutoff_1': cutoff_1_ls,
                                   'Cutoff_2': cutoff_2_ls,
                                   })
        summary_df.to_excel(f'../result/train_model/simplify/{model}_train_model_cv_{n_splits}_df.xlsx')
