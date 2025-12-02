import numpy as np
import pandas as pd
import gseapy as gp
import scanpy as sc
import pandas as pd
import os
import warnings
import argparse
from pathlib import Path


cli = argparse.ArgumentParser(description="Validate GSEA null scenarios.")
cli.add_argument("--null_type", choices=["null_overlap", "null_sparsity", "null_sample"], default="null_overlap", help="Select which null dataset to evaluate.")
args = cli.parse_args()


DATA_ROOT="/disk/tanbowen/scKnockPath/data/simulation"
RESULTS_ROOT="/disk/tanbowen/scKnockPath/results/simulation_exp"

null_type = args.null_type

if null_type == "null_overlap":
    data_subdir = "null_overlap_data"
    prefix= "overlap"
    results_subdir = "GSEA_null_overlap_results"
    param_values = range(0, 11, 2)
    
elif null_type == "null_sparsity":
    data_subdir = "null_sparsity_data"
    prefix= "sparsity"
    results_subdir = "GSEA_null_sparsity_results"
    param_values = range(6, 17, 2)
else:  # null_sample
    data_subdir = "null_sample_data"

    prefix= "sample"
    results_subdir = "GSEA_null_sample_results"
    param_values = range(600, 1801, 240)

seeds = range(46, 76)

def GSEA_test(filename,fdr=0.2,null_type=null_type):
    adata = sc.read_h5ad(filename)
    
    if null_type != "null_sample":
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()

    true_effect_pathways = adata.uns['effect_pathways']

    df = pd.DataFrame(adata.layers['lognorm'], columns=adata.var_names)
    geneset_dict = adata.uns['all_pathways']
    obs_y = 'cell_type'

    res = gp.gsea(
        data=df.T,
        gene_sets={pathway: genes.astype(str).tolist() for pathway, genes in geneset_dict.items()},
        cls=adata.obs[obs_y],
        permutation_num=1000,
        method='t_test',
        min_size=1,
        seed=37
    )
    fdr = 0.2
    # print(res.res2d)
    selected_pathways=res.res2d[(res.res2d['FDR q-val'] <= fdr)]['Term']
    
    return len(selected_pathways)


warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

n_selected_pathways_list = {}

for param in param_values:
    print(f'{null_type}={param}')
    path = f"{DATA_ROOT}/{data_subdir}/n_{prefix}_beta{param}"
    n_selected_pathways = []
    for seed in range(46, 76):
        print(f'seed={seed}')
        n_selected= GSEA_test(f"{path}/simu_scRNAseq_100pathways_{param}{prefix}_null_seed={seed}.h5ad")
        n_selected_pathways.append(n_selected)
        print(f'n_selected_pathways={n_selected}')
        print('--------------------------------------------------------------------')

    print(n_selected_pathways)

    n_selected_pathways_list[str(param)] = n_selected_pathways
    print('===============================================')

n_selected_pathways_df = pd.DataFrame(n_selected_pathways_list)
# Create a folder to store results
results_folder = f"{RESULTS_ROOT}/{results_subdir}"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("GSEA_overlap_results folder created!")

# Save results in the new folder
n_selected_pathways_df.to_csv(os.path.join(results_folder, f"GSEA_{null_type}_n_selected_values.csv"), index=False)

    