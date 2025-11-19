#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/27
# @Author  : Qiang Huang
# @File    :
'''
    shannon diversity for emb cluster，boxplot
'''
import os
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statannotations.Annotator import Annotator

if __name__ == '__main__':
    os.chdir('../../')
    datestamp = datetime.now()
    date = f'{datestamp.year}{datestamp.month:02d}{datestamp.day:02d}'
    print(date)
    # input
    shan_summary = pd.read_excel('../result/deepair/embedding/shannon_eval_250927.xlsx', index_col=0).T

    # boxplot
    Types = ['CDR3B', 'Vbeta', 'CDR3s', 'Vab']
    models = ['AF3', 'tfold-TCR', 'TCR-BERT']
    shan_sum_lon = shan_summary.melt()
    shan_sum_lon['Type'] = shan_sum_lon['mode'].map(lambda x: x.split('_')[1])
    shan_sum_lon['Type'] = pd.Categorical(shan_sum_lon['Type'],
                                           categories=Types,
                                           ordered=True)
    shan_sum_lon['model'] = shan_sum_lon['mode'].map(lambda x: x.split('_')[0])
    shan_sum_lon['model'] = pd.Categorical(shan_sum_lon['model'],
                                           categories=models,
                                           ordered=True)
    sig_pairs1 = [
        (('CDR3B', models[0]), ('CDR3B', models[2])),
        (('CDR3B', models[1]), ('CDR3B', models[2])),
        (('CDR3s', models[0]), ('CDR3s', models[2])),
        (('CDR3s', models[1]), ('CDR3s', models[2]))
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
    plt.xlabel('')
    plt.ylim(y_lims[0], 1.15 * y_lims[1])
    plt.xticks(ticks=[0, 1, 2, 3], labels=Types, fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel('Shannon Index', fontsize=20)
    plt.legend(loc='upper left', fontsize=24)
    plt.savefig(f'../result/deepair/embedding/leiden_cluster_peptide_shan_eval_{date}.pdf',
                dpi=300, bbox_inches='tight')