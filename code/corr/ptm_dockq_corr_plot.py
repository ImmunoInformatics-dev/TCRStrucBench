#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TCRstruc
# @Purpose : visualize the evaluation index
# @Time    : 2025/9/15
# @Author  : Qiang Huang
# @File    : 
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# global setting
plt.rcParams['axes.titlesize'] = 20 
plt.rcParams['axes.labelsize'] = 20  
plt.rcParams['xtick.labelsize'] = 16 
plt.rcParams['ytick.labelsize'] = 20 
plt.rcParams['legend.fontsize'] = 18   
plt.rcParams['legend.title_fontsize'] = 20  

def p_label(p_value):
    if p_value > 0.05:
        label = 'NoSig'
    elif p_value > 0.01:
        label = '*'
    elif p_value > 0.001:
        label = '**'
    elif p_value > 0.0001:
        label = '***'
    else:
        label = '****'
    return label

def corr_func(x, y, method):
    if method is None:
        _, pvalue1 = stats.shapiro(x)
        _, pvalue2 = stats.shapiro(y)
        if pvalue1 > 0.05 and pvalue2 > 0.05:
            method = 'Pearson'
            coef, pvalue = stats.pearsonr(x, y)
            label = p_label(pvalue)
        else:
            method = 'Spearman'
            coef, pvalue = stats.spearmanr(x, y)
            label = p_label(pvalue)
    elif method == 'Pearson':
        coef, pvalue = stats.pearsonr(x, y)
        label = p_label(pvalue)
    elif method == 'Spearman':
        coef, pvalue = stats.spearmanr(x, y)
        label = p_label(pvalue)
    else:
        raise ValueError('''method is wrong, try 'Pearson' or 'Spearman' ''')

    return coef, pvalue, label, method

def corr_stat(input_file1, input_file2, model_name, out_file,
              method=None):
    '''
    :param input_file1: plddt output file
    :param inputfile2: rmsd output file
    :param sheet_name: model names
    :param domain_ls:
    :param target:
    :param method: correlation methods: pearson or spearman
    :param plot: whether to plot the heatmap
    :param out_file: summary corr. table file
    :param fig_file: save the fig. of plot
    :return:
    '''
    corr_dict = {'Model': [], 'Coef': [], 'pvalue': [], 'label': [], 'method': []}

    ptm_df = pd.read_excel(input_file1, index_col=0)
    target_df = pd.read_excel(input_file2, index_col=0)
    target_df = target_df.loc[ptm_df["pdb_id"], :]
    for model in model_name:
        corr_dict['Model'].append(model)
        ptm_ls = ptm_df[model].values
        target_ls = target_df[model].values
        coef, pvalue, label, _method = corr_func(ptm_ls, target_ls, method=method)
        corr_dict['Coef'].append(coef)
        corr_dict['pvalue'].append(pvalue)
        corr_dict['label'].append(label)
        corr_dict['method'].append(_method)
    #
    corr_df = pd.DataFrame(corr_dict)
    corr_df.to_excel(out_file)

def scatter_plot(input_file1, input_file2, model_name, show=True, fig_file=None):
    '''
    :param input_file1: ptm output file
    :param inputfile2: rmsd output file
    :param sheet_name: model names
    :param method: correlation methods: pearson or spearman
    :param domain_ls: the domain to choose
    :param fig_file: save the fig. of plot
    :return:
    '''
    fig, axs = plt.subplots(2, 3, figsize=(36, 24))
    ptm_df = pd.read_excel(input_file1, index_col=0)
    target_df = pd.read_excel(input_file2, index_col=0)
    target_df = target_df.loc[ptm_df["pdb_id"], :]
    for i, ax in enumerate(axs.flat):
        if i < len(model_name):
            model = model_name[i]
            ptm_ls = ptm_df[model].values
            target_ls = target_df[model].values
            coef, p_value, label, method = corr_func(ptm_ls, target_ls, method=None)
            text_anno = f'Coef.:{coef:.2f} ({method}) \nP_value:{p_value:.3f} ({label})'
            sns.set_context(font_scale=1.5)
            sns.regplot(
                x=ptm_ls, y=target_ls, ax=ax,
                scatter_kws={'alpha': 0.6, 'color': 'teal'},
                line_kws={'color': 'darkorange'}
            )
            ax.text(
                x=0.05, y=0.85,
                s=text_anno,
                transform=ax.transAxes,
                fontsize=30
            )
            ax.set_title(model_name[i], fontsize=36)
            ax.tick_params(labelsize=28)
            ax.set_ylabel('Dockq', fontsize=28)
            ax.set_xlabel('Ranking Score', fontsize=28)
            ax.grid(True)
        else:
            ax.axis('off')
    if show:
        plt.show()
    if fig_file:
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    os.chdir('../../')
    corr_stat(input_file1='../result/ptm/ptm_score.xlsx',
              input_file2='../result/dockq/DockQ_Models_summary.xlsx',
              model_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
              out_file='../result/corr/corr_paired_ptm_dockq.xlsx',
              method=None)

    scatter_plot(input_file1='../result/ptm/ptm_score.xlsx',
                 input_file2='../result/dockq/DockQ_Models_summary.xlsx',
                 model_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
                 show=False,
                 fig_file='../result/corr/scatterplot_paired_ptm_dockq.pdf')

