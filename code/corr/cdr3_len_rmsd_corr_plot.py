#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TCRstruc
# @Purpose : visualize the evaluation index
# @Time    : 2025/9/22
# @Author  : Qiang Huang
# @File    : plot_250428.py
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

def corr_stat(cdr3_df, input_file, sheet_name, domain_ls, out_file,
              method=None):
    corr_dict = {'Model': [], 'Domain': [], 'Coef': [], 'pvalue': [], 'label': [], 'method': []}
    for sheet in sheet_name:
        cdr3_df_w = cdr3_df.pivot(index='PDB_ID', columns='Domain', values='CDR3_length')
        cdr3_df_w = cdr3_df_w.reset_index()
        target_df = pd.read_excel(input_file, sheet_name=sheet, index_col=1)
        target_df = target_df.loc[cdr3_df_w['PDB_ID'], :]
        for idx, domain in enumerate(domain_ls):
            corr_dict['Model'].append(sheet)
            corr_dict['Domain'].append(domain)
            cdr3_len = cdr3_df_w[domain].values
            target_ls = target_df[domain].values
            coef, pvalue, label, _method = corr_func(cdr3_len, target_ls, method=method)
            corr_dict['Coef'].append(coef)
            corr_dict['pvalue'].append(pvalue)
            corr_dict['label'].append(label)
            corr_dict['method'].append(_method)
    #
    corr_df = pd.DataFrame({'Model': list(corr_dict['Model']), 'Domain': list(corr_dict['Domain']), 'Coef': list(corr_dict['Coef'])})
    corr_df_w = corr_df.pivot(columns='Domain', index='Model', values='Coef')
    corr_df_w = corr_df_w.loc[sheet_name, domain_ls]

    p_df = pd.DataFrame(
        {'Model': list(corr_dict['Model']), 'Domain': list(corr_dict['Domain']), 'Annotation': corr_dict['pvalue']})
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
    return corr_df_w, label_df_w, mask_df_w

def scatter_plot(cdr3_df, input_file, sheet_name, domain_ls, target,
                 corr_df_w, mask_df_w, show=True, fig_file=None):
    palette = {"A-CDR3": "#66c2a5", "B-CDR3": "#fc8d62"}
    line_styles = {"A-CDR3": "--", "B-CDR3": "-"}
    fig, axs = plt.subplots(2, 3, figsize=(38, 20))
    for idx, ax in enumerate(axs.flat):
        if idx < len(sheet_name):
            sheet = sheet_name[idx]
            target_df = pd.read_excel(input_file, sheet_name=sheet)
            target_ls = target_df.loc[:, ['PDB_ID'] + domain_ls]
            target_lon = target_ls.melt(id_vars='PDB_ID', var_name='Domain', value_name='RMSD')
            target_final = pd.merge(target_lon, cdr3_df, on=['PDB_ID', 'Domain'])
            sns.set_context(font_scale=1.5)
            sns.scatterplot(
                data=target_final,
                x="CDR3_length",
                y="RMSD",
                hue="Domain",  
                style="Domain", 
                ax=ax,
                palette=palette,  
                s=100  
            )
            # 
            text_anno = []
            for d in domain_ls:
                target_subset = target_final[target_final["Domain"] == d]
                cdr3_len = target_subset["CDR3_length"].values
                target_ls = target_subset["RMSD"].values
                coef, p_value, label, method = corr_func(cdr3_len, target_ls, method=None)
                text_anno.append(f'{d} Coef.:{coef:.2f} ({label})')
                sns.regplot(
                    data=target_subset,
                    x="CDR3_length",
                    y="RMSD",
                    ax=ax,
                    scatter=False,
                    robust=True,
                    color=palette[d],
                    line_kws={
                        "linestyle": line_styles[d],
                        "alpha": 0.7
                    }
                )
            ax.text(  # A chain
                x=0.05, y=0.9,
                s=text_anno[0],
                transform=ax.transAxes,
                fontsize=26
            )
            ax.text(  # B chain
                x=0.05, y=0.8,
                s=text_anno[1],
                transform=ax.transAxes,
                fontsize=26
            )
            ax.set_title(sheet_name[idx], fontsize=36)
            ax.tick_params(labelsize=28)
            ax.set_ylabel(target, fontsize=28)
            ax.set_xlabel('CDR3 length', fontsize=28)
            ax.legend(loc='upper right', fontsize=28)
            ax.grid(True)
        else:
            sns.heatmap(
                corr_df_w,
                mask=mask_df_w,
                annot=True,
                vmax=0.7,
                vmin=0,
                fmt=".2f",
                annot_kws={"size": 30, "color": "black", "weight": "bold"},
                cmap='Greens', center=0.3,
                linewidths=0.5, linecolor="lightgray"
            )
            ax.set_title(f"Correlation between CDR3 length and {target}")

    if show:
        plt.show()
        plt.close()
    if fig_file:
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close()



if __name__ == '__main__':
    os.chdir('../../')
    tcr_df = pd.read_csv('../result/annotation/tcr_anno_df.csv')
    tcr_df['CDR3_length'] = tcr_df['CDR3'].map(lambda x: len(x))
    tcr_df['Domain'] = ['A-CDR3', 'B-CDR3'] * 51
    tcr_df['PDB_ID'] = tcr_df['Name'].map(lambda x: x.split('_')[0])
    cdr3_df = tcr_df.loc[:, ['PDB_ID', 'Domain', 'CDR3_length']]

    # plddt_rmsd - single
    corr_df_w, _, mask_df_w = corr_stat(cdr3_df=cdr3_df,
              input_file='../result/rmsd/rmsd_summary_single.xlsx',
              sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
              domain_ls=['A-CDR3', 'B-CDR3'],
              out_file='../result/corr/corr_single_cdr3_rmsd_0922.xlsx',
              method=None)

    scatter_plot(cdr3_df=cdr3_df,
                 input_file='../result/rmsd/rmsd_summary_single.xlsx',
                 sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
                 domain_ls=['A-CDR3', 'B-CDR3'],
                 target='RMSD',
                 corr_df_w=corr_df_w, mask_df_w=mask_df_w,
                 show=False,
                 fig_file='../result/corr/scatterplot_single_cdr3_rmsd.pdf')

    # plddt_rmsd - paired
    corr_df_w, _, label_df_w = corr_stat(cdr3_df=cdr3_df,
                                      input_file='../result/rmsd/rmsd_summary_paired.xlsx',
                                      sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
                                      domain_ls=['A-CDR3', 'B-CDR3'],
                                      out_file='../result/corr/corr_paired_cdr3_rmsd_0922.xlsx',
                                      method=None)

    scatter_plot(cdr3_df=cdr3_df,
                 input_file='../result/rmsd/rmsd_summary_paired.xlsx',
                 sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'],
                 domain_ls=['A-CDR3', 'B-CDR3'],
                 target='RMSD',
                 corr_df_w=corr_df_w, mask_df_w=mask_df_w,
                 show=False,
                 fig_file='../result/corr/scatterplot_paired_cdr3_rmsd.pdf')
