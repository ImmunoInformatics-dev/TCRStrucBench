#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project :
# @Purpose :
# @Time    : 2025/9/27
# @Author  : Qiang Huang
# @File    :
import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

def get_region_range(part_seq, whole_seq, intercept):
    '''
    '''
    start = whole_seq.find(part_seq[:5])  # first match
    end = start + len(part_seq)
    if end > len(whole_seq):
        end = len(whole_seq)
    start += intercept
    end += intercept
    return start, end

def extrac_emb(tcr_anno, emb_dir, file_suffix):
    for idx, (id, cdr3a, trav, cdr3b, trbv) in enumerate(zip(tcr_anno['ID'],
                                                tcr_anno['CDR3A'], tcr_anno['VAdomain'],
                                                tcr_anno['CDR3B'], tcr_anno['VBdomain'])):
        startA, endA = get_region_range(cdr3a, trav, 0)
        startB, endB = get_region_range(cdr3b, trbv, len(trav))
        file_name = f'{emb_dir}/id_{id}_{file_suffix}'
        emb = np.load(file_name)
        if 'npz' in file_suffix:
            emb = emb['single_embeddings']
        emb = emb.squeeze()
        if emb.shape[-2] > 60:  # full length of TRBV
            embA = emb[startA:endA, :]
            embB = emb[startB:endB, :]
            emb = np.concatenate([embA, embB], axis=0)
        emb_mean = np.mean(emb, axis=0)
        if idx == 0:
            emb_mtx = emb_mean
        else:
            emb_mtx = np.vstack((emb_mtx, emb_mean))
    return emb_mtx

#
if __name__ == "__main__":
    os.chdir('../../')
    tcr_anno_1 = pd.read_csv("../data/deepair/BRP_merging_final_250704.csv")
    tcr_anno_2 = pd.read_csv("../data/deepair/BRP_merging_final_250715.csv")
    tcr_anno = tcr_anno_2[[x in tcr_anno_1['ID'].values for x in tcr_anno_2['ID'].values]]
    model_names = ['tfold_tcr']
    suffix_names = ['sfea.npy']
    for m, s in zip(model_names, suffix_names):
        if m != 'AF3':
            model_emb = extrac_emb(tcr_anno=tcr_anno,
                                      emb_dir=f'../result/deepair/{m}_paired',
                                      file_suffix=f'{s}')
            np.save(f'../result/deepair/adata/{m}_cdr3s_emb_250927.npy', model_emb)
        else:
            model_emb = np.load(f'../result/deepair/{m}_cdr3s/{m}_cdr3s_emb_250927.npy')


        model_adata = ad.AnnData(X=model_emb, obs=tcr_anno)
        model_adata.obs_names = model_adata.obs['ID'].astype(str)
        model_adata.obs = model_adata.obs.drop(columns='ID')

        model_adata.write_h5ad(filename=f'../result/deepair/adata/{m}_adata_250927.h5ad')
        # PCA
        sc.pp.pca(model_adata, n_comps=50)
        sc.pp.neighbors(model_adata)
        sc.tl.umap(model_adata)
        sc.tl.tsne(model_adata)
        model_adata.write_h5ad(filename=f'../result/deepair/adata/{m}_adata_pca_250927.h5ad')






