#!/bin/bash
# CDR3B
mkdir -p ../result/AF3/foldseek_grid/foldseek_cdr3b_c_0.7_i_0.5_a_1/
mkdir -p ../result/tfold/foldseek_grid/foldseek_cdr3b_c_0.7_i_0.5_a_1/
foldseek easy-cluster ../result/AF3/deepair_cdr3b/ ../result/AF3/foldseek_grid/foldseek_cdr3b_c_0.7_i_0.5_a_1/ ../result/AF3/foldseek_af3_tmp -c 0.7 --alignment-type 1 --min-seq-id 0.5
foldseek easy-cluster ../result/tfold/deepair_cdr3b/ ../result/tfold/foldseek_grid/foldseek_cdr3b_c_0.7_i_0.5_a_1/ ../result/tfold/foldseek_tfold_tmp -c 0.7 --alignment-type 1 --min-seq-id 0.5
# CDR3S
mkdir -p ../result/AF3/foldseek_grid/foldseek_cdr3s_c_0.7_i_0.5_a_2/
mkdir -p ../result/tfold/foldseek_grid/foldseek_cdr3s_c_0.7_i_0.5_a_2/
foldseek easy-multimercluster ../result/AF3/deepair_paired_cdr3s/ ../result/AF3/foldseek_grid/foldseek_cdr3s_c_0.7_i_0.5_a_2/ ../result/AF3/foldseek_af3_tmp -c 0.7 --alignment-type 2 --min-seq-id 0.5 --multimer-tm-threshold 0.65 --chain-tm-threshold 0.5 --interface-lddt-threshold 0.65
foldseek easy-multimercluster ../result/tfold/deepair_paired_cdr3s/ ../result/foldseek_grid/tfold/foldseek_cdr3s_c_0.7_i_0.5_a_2/ ../result/tfold/foldseek_tfold_tmp -c 0.7 --alignment-type 2 --min-seq-id 0.5 --multimer-tm-threshold 0.65 --chain-tm-threshold 0.5 --interface-lddt-threshold 0.65

# Vbeta
mkdir -p ../result/AF3/foldseek_grid/foldseek_vbeta_c_0.7_i_0.9_a_2/
mkdir -p ../result/tfold/foldseek_grid/foldseek_vbeta_c_0.7_i_0.9_a_2/
foldseek easy-cluster ../result/AF3/deepair_vbeta/ ../result/AF3/foldseek_grid/foldseek_vbeta_c_0.7_i_0.9_a_2/ ../result/AF3/foldseek_af3_tmp -c 0.7 --alignment-type 2 --min-seq-id 0.9
foldseek easy-cluster ../result/tfold/deepair_vbeta/ ../result/tfold/foldseek_grid/foldseek_vbeta_c_0.7_i_0.9_a_2/ ../result/tfold/foldseek_tfold_tmp -c 0.7 --alignment-type 2 --min-seq-id 0.9

#Vab
mkdir -p ../result/AF3/foldseek_grid/foldseek_vab_c_0.7_i_0.6_a_1/
mkdir -p ../result/tfold/foldseek_grid/foldseek_vab_c_0.7_i_0.6_a_1/
foldseek easy-multimercluster ../result/AF3/deepair_vab/ ../result/AF3/foldseek_grid/foldseek_vab_c_0.7_i_0.6_a_1/ ../result/AF3/foldseek_af3_tmp -c 0.7 --alignment-type 1 --min-seq-id 0.6 --multimer-tm-threshold 0.65 --chain-tm-threshold 0.5 --interface-lddt-threshold 0.65
foldseek easy-multimercluster ../result/tfold/deepair_vab/ ../result/tfold/foldseek_grid/foldseek_vab_c_0.7_i_0.6_a_1/ ../result/tfold/foldseek_tfold_tmp -c 0.7 --alignment-type 1 --min-seq-id 0.6 --multimer-tm-threshold 0.65 --chain-tm-threshold 0.5 --interface-lddt-threshold 0.65