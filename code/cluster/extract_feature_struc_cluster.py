#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/25
# @Author  : Qiang Huang
# @File    :
'''
    extract IRS and DIE major cluster sequence for alignment
'''
import pandas as pd
import numpy as np
import os
import shutil
import logomaker
from weblogo import logo_formatter
from weblogo import logo
from Bio import AlignIO
import subprocess
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    tcr_anno = pd.read_csv("../data/deepair/BRP_merging_final_250704.csv")
    tcr_anno['ID'] = tcr_anno['ID'].astype(str)
    tcr_anno['Peptide'] = tcr_anno['Peptide'].map(peptide_dict)
    #
    id_peptide = pd.read_csv(f'../data/deepair/BRP_merging_final_250704.csv')
    id_peptide['Peptide'] = id_peptide['Peptide'].map(peptide_dict)
    cov, iden, aln = 0.7, 0.5, 1
    #
    model_names = ['AF3', 'tfold']
    pattern_ls = ['IRS', 'DIE']
    pep_ls = ['GIL_MP-Flu', 'YLQ_SARS-CoV-2']
    cluster_ls = ['id_113', 'id_22320']
    for model in model_names:
        # af3 - foldseek
        cluster_file = f'../result/deepair/foldseek_grid/foldseek_grid_{model}/foldseek_cdr3b_c_{cov}_i_{iden}_a_{aln}/_cluster.tsv'
        fs_cluster = pd.read_table(cluster_file, sep='\t', names=['cluster', 'ID'])
        fs_cluster['cluster'] = fs_cluster['cluster'].map(lambda x: x.rsplit('_',maxsplit=1)[0])
        fs_cluster['ID'] = fs_cluster['ID'].map(lambda x: int(x.split('_')[1]))
        peptide_cluster = pd.merge(id_peptide, fs_cluster, on='ID', how='left')
        for pep, clu in zip(pep_ls, cluster_ls):
            my_dir = f'../result/deepair/foldseek_{model}_clu/{clu}'
            os.makedirs(my_dir, exist_ok=True)
            peptide_cluster_extract = peptide_cluster[(peptide_cluster['cluster'] == clu) & (peptide_cluster['Peptide'] == pep)]
            cdr3b = list(peptide_cluster_extract['CDR3B'])
            my_id = list(peptide_cluster_extract['ID'])
            # msa
            with open(f'{my_dir}/{clu}.fasta', "w") as f:
                for idx, seq in zip(my_id, cdr3b):
                    f.write(f">seq_{idx}\n{seq}\n")
            input_file = f'D:/Z440/Postdoctoral/Research/TCRsructure/result/deepair/foldseek_{model}_clu/{clu}/{clu}.fasta'
            output_file = f'D:/Z440/Postdoctoral/Research/TCRsructure/result/deepair/foldseek_{model}_clu/{clu}/{clu}_aligned.fasta'
            # MUSCLE
            cmd = [
                "muscle",
                "-align", input_file,
                "-output", output_file
            ]
            subprocess.run(cmd,
                           check=True,
                           # stdout=subprocess.PIPE,
                           # stderr=subprocess.PIPE,
                           text=True,
                           encoding='utf-8')
            # align
            alignment = AlignIO.read(f'{my_dir}/{clu}_aligned.fasta', "fasta")
            aligned_sequences = [str(seq.seq) for seq in alignment]
            counts_df_x = logomaker.alignment_to_matrix(sequences=aligned_sequences, to_type='counts', characters_to_ignore='.-X')
            logo = logomaker.Logo(
                counts_df_x,
                color_scheme='chemistry',
                # font_name='Arial Rounded MT Bold',
                stack_order='big_on_top',  #
                figsize=(12, 4)
            )


            logo.ax.set_ylabel("Information (bits)")
            logo.ax.set_title(clu, fontsize=32)
            logo.ax.grid(False)
            logo.style_spines(visible=False)
            plt.tight_layout()
            plt.savefig(f"{my_dir}/{clu}_motif_logo.pdf", dpi=300)
            # structrue
            for idx in my_id:
                shutil.copy2(src=f'../result/{model}/deepair_cdr3b/id_{idx}_cdr3b.pdb', dst=f'{my_dir}/id_{idx}_cdr3b.pdb')
