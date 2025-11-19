#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : TCRstruc
# @Purpose : produce the code for different model
# @Time    : 2025/9/2
# @Author  : Qiang Huang
# @File    : batch_run_code.py
import os
import pandas as pd
import numpy as np
import json
# input
os.chdir('../')
tcr_anno = pd.read_csv("../result/annotation/tcr_anno_df.csv")
pdb_ids = tcr_anno['Name']
gpu_id = 'hgx'
# dict
tcr_seq_dic = {}
for kk, vv in zip(pdb_ids, tcr_anno['Vdomain']):
    tcr_seq_dic[kk] = vv

# wite fasta
os.makedirs('../data/TCR_Vdomain_test_single', exist_ok=True)
for pdb_id in tcr_seq_dic.keys():
    with open(f"../data/TCR_Vdomain_test_single/{pdb_id}.fasta", "w") as f:
        f.write(f'>{pdb_id} \n'
                f'{tcr_seq_dic[pdb_id]}')

## for TCRmodel2
with open(f"../code/tcrmodel2_code_test_single.lsf", "w") as f:
    f.write("#!/bin/bash\n")
    f.write("#BSUB -J tcrmodel2_test_single\n")
    f.write(f"#BSUB -q {gpu_id}\n")
    f.write("#BSUB -n 8\n")
    f.write(f"#BSUB -e tcrmodel2_test_single.err\n")
    f.write(f"#BSUB -o tcrmodel2_test_single.out\n")
    f.write('''#BSUB -R "span[ptile=8]"\n''')
    f.write('''#BSUB -gpu "num=1"\n'''
            '\n')
    for idx, pdb_id in enumerate(tcr_seq_dic.keys()):
        f.write('cd /work/USER/tcrmodel2 \n')
        f.write(f'mkdir -p /work/USER/TCRstruc/result/tcrmodel2/test_single/{pdb_id} \n')
        if idx % 2 == 0:  # odd is tra
            f.write(f'python run_tcrmodel2_ub_tcr.py \\\n'
                    f'  --job_id=single_{pdb_id} \\\n'
                    f'  --output_dir=/work/USER/TCRstruc/result/tcrmodel2/test_single/{pdb_id} \\\n'
                    f'  --tcra_seq={tcr_seq_dic[pdb_id]} \\\n'
                    '  --ori_db=/work/USER/alphafold-2.3.1/afdb \\\n'
                    '  --model_preset=monomer \\\n'
                    '  --models_to_relax=best \\\n'
                    '  --max_template_date=2023-01-31 \\\n')
        else:  # even is trb
            f.write(f'python run_tcrmodel2_ub_tcr.py \\\n'
                    f'  --job_id=single_{pdb_id} \\\n'
                    f'  --output_dir=/work/USER/TCRstruc/result/tcrmodel2/test_single/{pdb_id} \\\n'
                    f'  --tcrb_seq={tcr_seq_dic[pdb_id]} \\\n'
                    '  --ori_db=/work/USER/alphafold-2.3.1/afdb \\\n'
                    '  --model_preset=monomer \\\n'
                    '  --models_to_relax=best \\\n'
                    '  --max_template_date=2023-01-31 \\\n')
        f.write(f'cd /work/USER/TCRstruc/result/tcrmodel2/test_single/{pdb_id}/test_single_{pdb_id}/test_single_{pdb_id}\n'
                'ls | grep [0-9].pkl | xargs -i rm {}\n'
                'ls | grep unrelaxed | xargs -i rm {}\n'
                '\n')

## for AF2
with open(f"../code/af2_code_test_single.lsf", "w") as f:
    f.write("#!/bin/bash\n")
    f.write(f"#BSUB -J AF2_test_single\n")
    f.write(f"#BSUB -q {gpu_id}\n")
    f.write("#BSUB -n 8\n")
    f.write(f"#BSUB -e AF2_test_single.err\n")
    f.write(f"#BSUB -o AF2_test_single.out\n")
    f.write('''#BSUB -R "span[ptile=8]"\n''')
    f.write('''#BSUB -gpu "num=1"\n'''
            '\n')
    f.write(f'mkdir -p /work/USER/TCRstruc/result/AF2/test_single  \n')
    for pdb_id in tcr_seq_dic.keys():
        f.write("cd /work/USER/alphafold-2.3.1/ \n")
        f.write('bash run_alphafold.sh \\\n'
                '  -d /work/USER/alphafold-2.3.1/afdb/ \\\n'
                f'  -o /work/USER/TCRstruc/result/AF2/test_single \\\n'
                f'  -f /work/USER/TCRstruc/data/TCR_Vdomain_test_single/{pdb_id}.fasta \\\n'
                '  -t 2023-01-31 \\\n'
                '  -r best \\\n'
                '  -m monomer \\\n')
        f.write(f'cd /work/USER/TCRstruc/result/AF2/test_single/{pdb_id}\n'
                'ls | grep [0-9].pkl | xargs -i rm {}\n'
                'ls | grep unrelaxed | xargs -i rm {}\n'
                '\n')


## for AF3
os.makedirs(f'../data/AF3/test_single', exist_ok=True)
for pdb_id in tcr_seq_dic.keys():
    json_dict = {"name": f'{pdb_id}',
                 "sequences": [{"protein": {"id": "A", "sequence": f"{tcr_seq_dic[pdb_id]}"}}],
                 "modelSeeds": [42],
                 "dialect": "alphafold3",
                 "version": 1
                 }
    with open(f"../data/AF3/test_single/af3_test_single_{pdb_id}.json", "w", encoding="utf-8") as f:
        json.dump(json_dict, f, ensure_ascii=False, indent=4)
# run file
with open(f"../code/AF3/af3_code_test_single.lsf", "w") as f:
    f.write("#!/bin/bash\n")
    f.write(f"#BSUB -J AF3_test_single\n")
    f.write(f"#BSUB -q {gpu_id}\n")
    f.write("#BSUB -n 8\n")
    f.write(f"#BSUB -e AF3_test_single.err\n")
    f.write(f"#BSUB -o AF3_test_single.out\n")
    f.write('''#BSUB -R "span[ptile=8]"\n''')
    f.write('''#BSUB -gpu "num=1"\n'''
            '\n')
    f.write("module purge \n"
            "module load  cuda/12.6 \n"
            "source /share/apps/anaconda3/bin/activate /share/apps/anaconda3/envs/alphafold3-v3.0.1 \n"
            "\n"
            )
    f.write('''export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false" \n'''
            "export XLA_PYTHON_CLIENT_PREALLOCATE=true \n"
            "export XLA_CLIENT_MEM_FRACTION=0.95 \n"
            "\n"
            )

    f.write(f'mkdir -p /work/USER/TCRstruc/result/AF3/test_single \n'
            "\n")

    for pdb_id in tcr_seq_dic.keys():
        f.write('python /share/apps/alphafold3-v3.0.1/run_alphafold.py \\\n'
                f'  --json_path=/work/USER/TCRstruc/data/AF3/test_single/af3_test_single_{pdb_id}.json \\\n'
                '  --model_dir=/work/USER/AF3 \\\n'
                f'  --db_dir=/share/apps/alphafold3-data \\\n'
                '  --max_template_date=2023-01-31 \\\n'
                '  --save_embeddings=True \\\n'
                f'  --output_dir=/work/USER/TCRstruc/result/AF3/test_single > $LSB_JOBID.log 2>&1 \n'
                '\n'
                )

## for esmfold2
with open(f'esmfold2_batch_predict_test_single.py', 'w') as f:
    f.write('import torch, esm, os \n'
            'import numpy as np \n'
            '\n'
            'model = esm.pretrained.esmfold_v1() \n'
            'model = model.eval().cuda() \n'
            '''os.makedirs('/work/USER/TCRstruc/result/ESMFold2/test_single', exist_ok=True) \n'''
            '''os.makedirs('/work/USER/TCRstruc/result/ESMFold2/test_single_ptm', exist_ok=True) \n'''
            '\n')

    for pdb_id in tcr_seq_dic.keys():
        f.write(f'''pdb_file = "/work/USER/TCRstruc/result/ESMFold2/test_single/{pdb_id}.pdb" \n'''
                f'''ptm_file = "/work/USER/TCRstruc/result/ESMFold2/test_single_ptm/{pdb_id}.npy" \n'''
                f'''fasta_seq = "{tcr_seq_dic[pdb_id]}" \n'''
                'with torch.no_grad(): \n'
                '\tpdbs, _, ptm, _ = model.infer_pdb(fasta_seq) \n'
                f'''with open(pdb_file, "w") as f: \n'''
                '\tf.write(pdbs) \n'
                'np.save(ptm_file, ptm.cpu().numpy()) \n'
                '\n')

## for tfold-TCR
os.makedirs('../data/tfold/test_single', exist_ok=True)
for pdb_id in tcr_seq_dic.keys():
    json_dict = {
                    "name": f'{pdb_id}',
                     "chains": [
                         {
                             "id": "A",
                             "sequence": f"{tcr_seq_dic[pdb_id]}"
                        }
                     ]
                 }
    with open(f"../data/tfold/test_single/tfold_test_single_{pdb_id}.json", "w", encoding="utf-8") as f:
        json.dump([json_dict], f, ensure_ascii=False, indent=4)


with open("../code/tfold-TCR/tfold_code_test_single.sh", "w") as f:
    f.write("#!/bin/bash\n"
            "\n")
    f.write(f'mkdir -p /mnt/zhangzheng_group/huangq-54/TCRstrc/result/tfold/test_single \n'
            "\n")
    for pdb_id in tcr_seq_dic.keys():
        f.write('python /mnt/zhangzheng_group/huangq-54/tfold/projects/tfold_tcr/predict.py \\\n'
                f'  --json=/mnt/zhangzheng_group/huangq-54/TCRstrc/data/tfold/test_single/tfold_test_single_{pdb_id}.json \\\n'
                f'  --output=/mnt/zhangzheng_group/huangq-54/TCRstrc/result/tfold/test_single/ \\\n'
                '  --model_version=TCR \\\n'
                f'  --seed=42 \n'
                '\n'
                )
