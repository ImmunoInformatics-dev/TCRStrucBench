#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/25
# @Author  : Qiang Huang
# @File    :
'''
    shannon diversity for struc cluster, boxplot
'''
import os
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statannotations.Annotator import Annotator

def shannon_ind(x):
    x_prop = x/sum(x)
    shan = x_prop * np.log(x_prop)
    shan_index = - shan.sum()
    return shan_index
#
peptide_dict = {
    'ELAGIGILTV': 'ELA_MART-1_Can',
    'GILGFVFTL': 'GIL_MP-Flu',
    'GLCTLVAML': 'GLC_BMLF1_EBV',
    'KLGGALQAK': 'KLG_IE-1_CMV',
    'AVFDRKSDAK': 'AVF_EBNA-3B_EBV',
    'IVTDFSVIK': 'IVT_EBNA-3B_EBV',
    'RAKFKQLL': 'RAK_BZLF1_EBV',
    'LTDEMIAQY': 'LTD_SARS-CoV-2',
    'TTDPSFLGRY': 'TTD_SARS-CoV-2',
    'YLQPRTFLL': 'YLQ_SARS-CoV-2',
}

if __name__ == '__main__':
    os.chdir('../../')
    datestamp = datetime.now()
    date = f'{str(datestamp.year)[2:]}{datestamp.month:02d}{datestamp.day:02d}'
    print(date)
    types = ['cdr3b', 'vbeta', 'cdr3s', 'vab']
    model_names = ['AF3', 'tfold']

    sh_dict = {}

    for tp in types:
        if tp in ['cdr3s', 'vab']:
            id_peptide = pd.read_csv(f'../data/deepair/BRP_merging_final_250715.csv')
        else:
            id_peptide = pd.read_csv(f'../data/deepair/BRP_merging_final_250704.csv')
        id_peptide['Peptide'] = id_peptide['Peptide'].map(peptide_dict)

        if tp == 'cdr3b':
            cov, iden, aln = 0.7, 0.5, 1
        elif tp == 'cdr3s':
            cov, iden, aln = 0.7, 0.5, 2
        elif tp == 'vbeta':
            cov, iden, aln = 0.7, 0.9, 2
        else:
            cov, iden, aln = 0.7, 0.6, 1

        for model in model_names:
            cluster_file = f'../result/deepair/foldseek_grid/foldseek_grid_{model}/foldseek_{tp}_c_{cov}_i_{iden}_a_{aln}/_cluster.tsv'
            fs_cluster = pd.read_table(cluster_file, sep='\t', names=['cluster', 'ID'])
            fs_cluster['ID'] = fs_cluster['ID'].map(lambda x: int(x.split('_')[1]))
            cluster_count = fs_cluster['cluster'].value_counts()
            minor_ls = cluster_count.index.values[cluster_count < 5]
            fs_cluster['cluster'] = ['minor_cluster' if x in minor_ls else x for x in fs_cluster['cluster']]

            peptide_cluster = pd.merge(id_peptide, fs_cluster, on='ID', how='left')
            peptide_cluster_filter = peptide_cluster.dropna(subset=['cluster'])
            #
            peptides = peptide_cluster_filter['Peptide']
            fd_cluster = peptide_cluster_filter['cluster']
            cm = pd.crosstab(fd_cluster, peptides)
            # shannon
            sh = cm.apply(axis=0, func=shannon_ind)
            sh_dict[f'{tp}_{model}'] = sh
    # add Gliph2 result
    cluster_file = '../result/deepair/tcr_anno_gliph2_250722.xlsx'
    gl_cluster = pd.read_excel(cluster_file)
    gl_cluster = gl_cluster.rename(columns={'pattern': 'cluster'})
    gl_cluster['cluster'] = ['minor_cluster' if x in ['single', 'no_pattern'] else x for x in gl_cluster['cluster']]
    cluster_count = gl_cluster['cluster'].value_counts()
    minor_ls = cluster_count.index.values[cluster_count < 5]
    gl_cluster['cluster'] = ['minor_cluster' if x in minor_ls else x for x in gl_cluster['cluster']]

    peptides = gl_cluster['Peptide']
    gl_cluster_ls = gl_cluster['cluster']
    cm = pd.crosstab(gl_cluster_ls, peptides)

    sh = cm.apply(axis=0, func=shannon_ind)
    sh_dict['cdr3b_Gliph2'] = sh

    # merge
    shan_summary = pd.DataFrame(sh_dict, index=cm.columns)
    shan_summary.to_excel(f'../result/deepair/foldseek_grid/foldseek_cluster_peptide_shan_eval_{date}.xlsx')

    Types = ['cdr3b', 'vbeta', 'cdr3s', 'vab']
    models = ['AF3', 'tfold', 'Gliph2']
    # boxplot
    shan_sum_lon = shan_summary.melt()
    shan_sum_lon['Type'] = shan_sum_lon['variable'].map(lambda x: x.split('_')[0])
    shan_sum_lon['Type'] = pd.Categorical(shan_sum_lon['Type'],
                                           categories=Types,
                                           ordered=True)
    shan_sum_lon['model'] = shan_sum_lon['variable'].map(lambda x: x.split('_')[1])

    # plot
    sig_pairs1 = [
        (('cdr3b', models[0]), ('cdr3b', models[2])),
        (('cdr3b', models[1]), ('cdr3b', models[2]))
    ]
    sig_pairs2 = [((tp, models[0]), (tp, models[1])) for tp in Types]
    sig_pairs = sig_pairs1 + sig_pairs2

    plt.figure(figsize=(8, 6))
    violin = sns.violinplot(data=shan_sum_lon,
                   x='Type',
                   y='value',
                   hue='model',
                   linewidth=0.5
                )
    annotator = Annotator(ax=violin, pairs=sig_pairs, data=shan_sum_lon, x='Type',
                          y='value', hue='model', order=Types)
    annotator.configure(test='Wilcoxon', text_format='star', loc='outside',
                        line_height=0.03, line_width=1, hide_non_significant=True)
    annotator.apply_and_annotate()

    y_lims = violin.get_ylim()
    plt.ylim(y_lims[0], 1.3 * y_lims[1])
    plt.xlabel('')
    plt.xticks(ticks=[0, 1, 2, 3], labels=Types, fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel('Shannon Index', fontsize=20)
    plt.legend(loc='lower right', fontsize=24)
    plt.savefig(f'../result/deepair/foldseek_grid/foldseek_cluster_peptide_shan_eval_{date}.pdf',
                dpi=300, bbox_inches='tight')