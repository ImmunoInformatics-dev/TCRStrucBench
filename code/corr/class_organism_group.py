#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/22
# @Author  : Qiang Huang
# @File    :
import itertools
import os

import pandas as pd
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
import json

def data_extract(pair_rmsd_file, single_rmsd_file, dockq_file, sheet_name):
    dockq = pd.read_excel(dockq_file)
    dockq = dockq.rename(columns={'Unnamed: 0': 'PDB_ID'})
    dockq_lon = dockq.melt(id_vars='PDB_ID', var_name='models')
    dockq_lon['domain'] = 'DockQ'
    Vab_p_df = pd.DataFrame()
    Vab_s_df = pd.DataFrame()
    for idx, s in enumerate(sheet_name):
        Vab_p_dfx = pd.read_excel(pair_rmsd_file, sheet_name=s, index_col=0).loc[:, ['PDB_ID', 'Valpha', 'Vbeta']]
        Vab_p_dfx['models'] = s
        Vab_s_dfx = pd.read_excel(single_rmsd_file, sheet_name=s, index_col=0).loc[:, ['PDB_ID', 'Valpha', 'Vbeta']]
        Vab_s_dfx['models'] = s
        if idx == 0:
            Vab_p_df = Vab_p_dfx
            Vab_s_df = Vab_s_dfx
        else:
            Vab_p_df = pd.concat([Vab_p_df, Vab_p_dfx], axis=0)
            Vab_s_df = pd.concat([Vab_s_df, Vab_s_dfx], axis=0)

    Vab_p_df_lon = Vab_p_df.melt(id_vars=['PDB_ID', 'models'], var_name='domain')
    Vab_p_df_lon['domain'] = Vab_p_df_lon['domain'] + '_paired'
    Vab_s_df_lon = Vab_s_df.melt(id_vars=['PDB_ID', 'models'], var_name='domain')
    Vab_s_df_lon['domain'] = Vab_s_df_lon['domain'] + '_single'
    final_out = pd.concat([Vab_p_df_lon, Vab_s_df_lon, dockq_lon], axis=0)
    return final_out

def box_plot_window(data,  x_lab, ylab, y, row_num, col_num, figsize, models,
                    categories, index_name, orders, sig_pairs=None,
                    outlier=True, jitter=False, fig_path=None):
    '''
    :param data:  dataframe including id, sample number, domain, rmsd
    :param group_col: for window group, domain col
    :param row_num: the figure rows
    :param col_num: the figure columns
    :param figsize: depend on your row and columns, eg. (width, height)
    :param groups: for each window names, domain names
    :return:
    '''
    fig, axes = plt.subplots(nrows=row_num, ncols=col_num, figsize=figsize)
    plt.subplots_adjust(wspace=0.25, hspace=0.3)  # 调整子图间距

    palette = sns.color_palette("Set2", n_colors=len(orders))

    for i, n in enumerate(index_name):
        for j, m in enumerate(models):
            window_data = data[(data['domain'] == n) & (data['models'] == m)]
            window_data['group'] = pd.Categorical(window_data['group'],
                                                  categories=orders,
                                                  ordered=True)

            # 使用seaborn绘制箱线图
            sns.boxplot(
                x=x_lab,
                y=y,
                data=window_data,
                ax=axes[i, j],
                palette=palette,
                width=0.6,
                showfliers=outlier,  
                notch=True,  
                showmeans=False,  
                meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"}
            )

            if jitter:
                sns.stripplot(
                    x=x_lab,
                    y=y,
                    data=window_data,
                    ax=axes[i, j],
                    color='black',
                    alpha=0.3,
                    jitter=0.2,
                    size=3
                )
            # group compare
            if sig_pairs:
                annotator = Annotator(axes[i, j], sig_pairs, data=window_data, x=x_lab, y=y, order=orders)
                annotator.configure(test='Levene', comparisons_correction='Bonferroni',
                                    text_format='star', line_height=0.03, line_width=1, hide_non_significant=True)
                annotator.apply_and_annotate()
            #
            y_lims = axes[i, j].get_ylim()
            axes[i, j].set_ylim(y_lims[0], 1.15 * y_lims[1])
            axes[i, j].set_title(f'{n}_{m}', fontsize=16, fontweight='bold')
            axes[i, j].set_xticklabels(categories, rotation=45, fontsize=14, ha='right')
            axes[i, j].set_xlabel('', fontsize=14)
            axes[i, j].set_ylabel(ylab[i], fontsize=14)
            axes[i, j].tick_params(axis='both', labelsize=14, which='major', direction='in', length=3, width=1., bottom=False)

            # 
            axes[i, j].grid(axis='y', linestyle='--', alpha=0.7)


    handles = []
    for color in palette:
        handles.append(plt.Rectangle((0, 0), 1, 1, color=color))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    if fig_path is not None:
        plt.savefig(fig_path, dpi=300)
        plt.close()

if __name__ == '__main__':
    os.chdir('../../')
    group_info = pd.read_excel('../data/TCR_Group_Info.xlsx')
    group_info = group_info.rename(columns={'Unnamed: 0': 'PDB_ID'})
    group_info = group_info.dropna(axis='rows', how='any')
    extract_df = data_extract(pair_rmsd_file='../result/rmsd/rmsd_summary_paired.xlsx',
                              single_rmsd_file='../result/rmsd/rmsd_summary_single.xlsx',
                              dockq_file='../result/dockq/DockQ_Models_summary.xlsx',
                              sheet_name=['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR'])
    pdb_ids = list(group_info['PDB_ID'])
    extract_df_filter = extract_df[[x in pdb_ids for x in extract_df['PDB_ID']]]
    final_df = pd.merge(extract_df_filter, group_info, on='PDB_ID')
    final_df['group'] = final_df['Ligand'].astype(str) + '_' + final_df['Organism'].astype(str)

    # Compare
    index_name = ['DockQ', 'Valpha_paired', 'Vbeta_paired', 'Valpha_single', 'Vbeta_single']
    model_names = ['AF2', 'TCRmodel2', 'AF3', 'ESMFoldv1', 'tfold-TCR']
    groups = ['ClassI_Homo', 'ClassI_Mus', 'ClassⅡ_Homo', 'ClassⅡ_Mus']
    # dif pairs
    sig_pairs = list(itertools.combinations(groups, 2))
    # plot
    box_plot_window(final_df, x_lab='group', models=model_names,
                    ylab=['DockQ', 'RMSD', 'RMSD', 'RMSD', 'RMSD'], y='value',
                    row_num=5, col_num=5, figsize=(24, 20),
                    categories=groups, index_name=index_name,
                    sig_pairs=sig_pairs, orders=groups,
                    jitter=True,
                    fig_path='../result/rmsd/HLA_Organism_group.pdf')