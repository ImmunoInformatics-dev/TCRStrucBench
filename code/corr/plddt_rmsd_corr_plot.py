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

# global
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

def corr_stat(input_file1, input_file2, sheet_name, domain_ls, target, out_file,
              method=None, plot=False, show=True, fig_file=None):
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
    corr_dict = {'Model': [], 'Domain': [], 'Coef': [], 'pvalue': [], 'label': [], 'method': []}

    for sheet in sheet_name:
        plddt_df = pd.read_excel(input_file1, sheet_name=sheet)
        target_df = pd.read_excel(input_file2, sheet_name=sheet)
        for domain in domain_ls:
            corr_dict['Model'].append(sheet)
            corr_dict['Domain'].append(domain)
            plddt_ls = plddt_df[domain].values
            target_ls = target_df[domain].values
            coef, pvalue, label, _method = corr_func(plddt_ls, target_ls, method=method)
            corr_dict['Coef'].append(coef)
            corr_dict['pvalue'].append(pvalue)
            corr_dict['label'].append(label)
            corr_dict['method'].append(_method)
    #
    corr_df = pd.DataFrame({'Model': list(corr_dict['Model']), 'Domain': list(corr_dict['Domain']), 'Coef': list(corr_dict['Coef'])})
    corr_df_w = corr_df.pivot(columns='Domain', index='Model', values='Coef')
    corr_df_w = corr_df_w.loc[sheet_name, domain_ls]

    p_df = pd.DataFrame({'Model': list(corr_dict['Model']), 'Domain': list(corr_dict['Domain']), 'Annotation': corr_dict['pvalue']})
    p_df_w = p_df.pivot(columns='Domain', index='Model', values='Annotation')
    p_df_w = p_df_w.loc[sheet_name, domain_ls]
    mask_df_w = p_df_w > 0.05

    label_df = pd.DataFrame({'Model': list(corr_dict['Model']), 'Domain': list(corr_dict['Domain']), 'Annotation': corr_dict['label']})
    label_df_w = label_df.pivot(columns='Domain', index='Model', values='Annotation')
    label_df_w = label_df_w.loc[sheet_name, domain_ls]

    with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
        corr_df_w.to_excel(writer, sheet_name='Coef.')
        p_df_w.to_excel(writer, sheet_name='P_value')
        label_df_w.to_excel(writer, sheet_name='Sig.')
    #
    if plot:
        plt.figure(figsize=(18, 3))
        sns.heatmap(
            corr_df_w,
            mask=mask_df_w,
            # annot=label_df_w,  
            annot=True,  
            fmt=".2f",  
            annot_kws={"size": 16, "color": "white", "weight": "bold"},
            # cmap="jet",
            vmin=-0.7, vmax=0,
            cmap="Greens_r",
            linewidths=0.5, linecolor="lightgray"
        )
        plt.title(f"Correlation Heatmap of pLDDT and {target}")
        if show:
            plt.show()
        if fig_file:
            plt.savefig(fig_file, dpi=300, bbox_inches='tight')

def scatter_plot(input_file1, input_file2, sheet_name, domain_ls, target, show=True, fig_file=None):
    '''
    :param input_file1: plddt output file
    :param inputfile2: rmsd output file
    :param sheet_name: model names
    :param method: correlation methods: pearson or spearman
    :param domain_ls: the domain to choose
    :param fig_file: save the fig. of plot
    :return:
    '''
    fig, axs = plt.subplots(len(sheet_name), len(domain_ls), figsize=(len(domain_ls)*10, len(sheet_name)*11))
    for i, sheet in enumerate(sheet_name):
        plddt_df = pd.read_excel(input_file1, sheet_name=sheet)
        target_df = pd.read_excel(input_file2, sheet_name=sheet)

        for j, domain in enumerate(domain_ls):
            plddt_ls = plddt_df[domain].values
            target_ls = target_df[domain].values
            coef, p_value, label, method = corr_func(plddt_ls, target_ls, method=None)
            text_anno = f'Coef.:{coef:.2f} ({method}) \nP_value:{p_value:.3f} ({label})'
            sns.set_context(font_scale=1.5)
            sns.regplot(
                x=plddt_ls, y=target_ls, ax=axs[i, j],
                scatter_kws={'alpha': 0.6, 'color': 'teal'},
                line_kws={'color': 'darkorange'}
            )
            axs[i, j].text(
                x=0.05, y=0.85,
                s=text_anno,
                transform=axs[i, j].transAxes,
                fontsize=30,
                # bbox={"facecolor": "white", "alpha": 0.8, "pad": 5}
            )
            axs[i, j].set_title(sheet + '_' + domain, fontsize=36)
            axs[i, j].tick_params(labelsize=28)
            axs[i, j].set_ylabel('rmsd', fontsize=28)
            axs[i, j].set_xlabel('plddt', fontsize=28)
            axs[i, j].grid(True)

    if show:
        plt.show()
    if fig_file:
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    os.chdir('../../')
    # plddt_rmsd - paired
    corr_stat(input_file1='../result/plddt/plddt_paired_summary.xlsx',
              input_file2='../result/rmsd/rmsd_summary_paired.xlsx',
              sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
              domain_ls=['Valpha', 'A-FWs', 'A-CDR1', 'A-CDR2', 'A-CDR3', 'Vbeta', 'B-FWs', 'B-CDR1', 'B-CDR2','B-CDR3'],
              target='RMSD',
              out_file='../result/corr/corr_paired_plddt_rmsd.xlsx',
              method=None,
              plot=True,
              show=False,
              fig_file='../result/corr/corr_paired_plddt_rmsd.pdf')

    scatter_plot(input_file1='../result/plddt/plddt_paired_summary.xlsx',
                 input_file2='../result/rmsd/rmsd_summary_paired.xlsx',
                 sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
                 domain_ls=['Valpha', 'A-FWs', 'A-CDR1', 'A-CDR2', 'A-CDR3', 'Vbeta', 'B-FWs', 'B-CDR1', 'B-CDR2','B-CDR3'],
                 target='RMSD',
                 show=False,
                 fig_file='../result/corr/scatterplot_paired_plddt_rmsd.pdf')
    # plddt_rmsd - single
    corr_stat(input_file1='../result/plddt/plddt_single_summary.xlsx',
              input_file2='../result/rmsd/rmsd_summary_single.xlsx',
              sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
              domain_ls=['Valpha', 'A-FWs', 'A-CDR1', 'A-CDR2', 'A-CDR3', 'Vbeta', 'B-FWs',  'B-CDR1', 'B-CDR2', 'B-CDR3'],
              target='RMSD',
              out_file='../result/corr/corr_single_plddt_rmsd.xlsx',
              method=None,
              plot=True,
              show=False,
              fig_file='../result/corr/corr_single_plddt_rmsd.pdf')

    scatter_plot(input_file1='../result/plddt/plddt_single_summary.xlsx',
                 input_file2='../result/rmsd/rmsd_summary_single.xlsx',
                 sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
                 domain_ls=['Valpha', 'A-FWs', 'A-CDR1', 'A-CDR2', 'A-CDR3', 'Vbeta', 'B-FWs',  'B-CDR1', 'B-CDR2', 'B-CDR3'],
                 target='RMSD',
                 show=False,
                 fig_file='../result/corr/scatterplot_single_plddt_rmsd.pdf')
