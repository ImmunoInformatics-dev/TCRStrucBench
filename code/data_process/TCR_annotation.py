#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TCRstruc
# @Purpose : annotate the TCR aa sequence with anarci
# @Time    : 2025/9/02
# @Author  : Qiang Huang
# @File    : TCR_aa_annotation.py
import os
from Bio import SeqIO
import subprocess
from anarci import anarci
import pandas as pd

def fasta_extract(fasta_file):
    '''
    :param fasta_file: fasta file
    :return: dictionary with id and sequence
    '''
    pdb_id = fasta_file[-10:-6]
    sequence_dic = {}
    for chain in SeqIO.parse(fasta_file, 'fasta'):
        chain_name = '_'.join([pdb_id, chain.description.split('|')[2].split(' ')[0]])  # ensure your sequence id format
        sequence_dic[chain_name] = str(chain.seq)
    return sequence_dic

def parse_anarci(in_seq):
    '''
    numbering the CDR3 and FR4 use Aho, and the others use IMGT
    :param in_seq: aa sequence
    :return:
    '''
    try:
        command_aho = f"ANARCI --scheme aho -i {in_seq}"
        command_imgt = f"ANARCI --scheme imgt -i {in_seq}"
        output_aho = subprocess.check_output(command_aho, shell=True, stderr=subprocess.STDOUT)
        output_imgt = subprocess.check_output(command_imgt, shell=True, stderr=subprocess.STDOUT)
        output_aho = output_aho.decode("utf-8")  # decode bytes to string
        output_imgt = output_imgt.decode("utf-8")  # decode bytes to string
    except subprocess.CalledProcessError as e:
        print(f"ANARCI failed for {in_seq} with error: {e.output.decode('utf-8')}")
        return "NA", "NA"

    tcr_dict = {'FW1': '', 'CDR1': '', 'FW2': '', 'CDR2': '', 'FW3': '', 'CDR3': '', 'FW4': '', 'Vdomain': ''}
    # imgt
    for i in output_imgt.split("\n"):
        if i and i[0] != "#" and i[0] != "/" and i[0] != " ":
            # _, num, res = i.rstrip().split()
            _, num, res = i.rstrip().split()[0:3]  # sometimes there are two AA at one location for imgt
            num = int(num)
            if res != "-":
                tcr_dict['Vdomain'] += res
                if num <= 26:
                    tcr_dict['FW1'] += res
                elif num >= 27 and num <= 38:
                    tcr_dict['CDR1'] += res
                elif num >= 39 and num <= 55:
                    tcr_dict['FW2'] += res
                elif num >= 56 and num <= 65:
                    tcr_dict['CDR2'] += res
                elif num >= 66 and num <= 104:
                    tcr_dict['FW3'] += res
    # Aho
    for i in output_aho.split("\n"):
        if i and i[0] != "#" and i[0] != "/" and i[0] != " ":
            _, num, res = i.rstrip().split()
            num = int(num)
            if res != "-":
                # tcr_dict['Vdomain'] += res
                if num >= 106 and num <= 139:
                    tcr_dict['CDR3'] += res
                elif num >= 140:
                    tcr_dict['FW4'] += res
    return tcr_dict


def get_germlines(seq:str):
    '''
    Get the VJ germlines from TCRa or TCRb sequences
    '''
    input_seq = [('seq', seq)]
    try:
        results = anarci(input_seq, scheme="aho", output=False, assign_germline=True)
        v_gene = results[1][0][0]['germlines']['v_gene'][0][1]
        j_gene = results[1][0][0]['germlines']['j_gene'][0][1]
    except:
        v_gene = 'NA'
        j_gene = 'NA'
    return v_gene, j_gene


def parse_tcr(in_seq):
    tcr_dict = parse_anarci(in_seq)
    v_gene, j_gene = get_germlines(in_seq)
    tcr_dict.update({'V_gene': v_gene})
    tcr_dict.update({'J_gene': j_gene})
    return tcr_dict

def update_vdomain(vdomain, wholeseq, match_len=5):
    '''
    not all vdomain are same with the origin tcr seq after anarci numbering
    :param vdomain: the calculated vdomain from anarci
    :param match_len: the start and end number of AA for matching
    :param wholeseq: the corresponding original tcr sequence
    :return: updated vdomain sequence
    '''
    pre_mer = vdomain[:match_len]
    end_mer = vdomain[-match_len:]
    start = wholeseq.find(pre_mer)  # first match
    end = wholeseq.find(end_mer) + len(end_mer)
    return wholeseq[start:end]


if __name__ == '__main__':
    fasta_files = os.listdir('../data/modelling_fasta')
    sequence_dict = {}
    for i in fasta_files:
        pdb_id = i.split('.')[0]
        # print(pdb_id)
        fasta_file = os.path.join('../data/modelling_fasta', i)
        sequence_dict.update(fasta_extract(fasta_file))

    tcr_dicts = {'Name': [], 'V_gene': [], 'J_gene': [], 'FW1': [], 'CDR1': [], 'FW2': [], 'CDR2': [], 'FW3': [], 'CDR3': [], 'FW4': [], 'Vdomain': []}

    for k, s in sequence_dict.items():
        tcr_dicts['Name'].append(k)
        tcr_anno = parse_tcr(s)
        # evaluate the agreement and update the vdomain
        if tcr_anno['Vdomain'] not in s:
            tcr_anno['Vdomain'] = update_vdomain(tcr_anno['Vdomain'], s, 5)
        for kk, ss in tcr_anno.items():
            tcr_dicts[kk].append(ss)

    tcr_anno_df = pd.DataFrame(tcr_dicts)
    tcr_anno_df.to_csv('../result/tcr_anno_df.csv', index=False)

