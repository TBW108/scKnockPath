import numpy as np
import pandas as pd
import gseapy as gp
import scanpy as sc
import pandas as pd
import os
import warnings


def GSEA_test(filename):
    adata = sc.read_h5ad(filename)
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
    fdrs = [0.2]
    # print(res.res2d)
    for fdr in fdrs:
        result = np.array(res.res2d[(res.res2d['FDR q-val'] <= fdr)]['Term'])
        real_fdr = np.logical_not(np.isin(result, true_effect_pathways)).sum() / len(result)
        real_power = np.isin(true_effect_pathways, result).sum() / len(true_effect_pathways)
        print(f'real_fdr={real_fdr}')
        print(f'real_power={real_power}')
        print(f'There are {len(result)} pathways selected with FDR<={fdr}')
    
    return real_fdr, real_power



warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

fdr_list = {}
power_list = {}

for n_overlap in range(0,11,2):
    print(f'n_overlap={n_overlap}')
    path = f"/disk/tanbowen/scKnockPath/data/simulation/overlap_data/n_overlap_beta{n_overlap}"
    fdr_values = []
    power_values = []
    
    for seed in range(46, 76):
        print(f'seed={seed}')
        fdr, power = GSEA_test(f"{path}/simu_scRNAseq_100pathways_{n_overlap}overlaplog_seed={seed}.h5ad")
        fdr_values.append(fdr)
        power_values.append(power)
        print(f'real_fdr={fdr}, real_power={power}')
        print('--------------------------------------------------------------------')
            
    print(fdr_values)
    print(power_values)
    
    fdr_list[str(n_overlap)] = fdr_values
    power_list[str(n_overlap)] = power_values
    print('===============================================')

fdr_df = pd.DataFrame(fdr_list)
power_df = pd.DataFrame(power_list)
# Create a folder to store results
results_folder = "/disk/tanbowen/scKnockPath/results/simulation_exp/GSEA_overlap_results"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("GSEA_overlap_results folder created!")

# Save results in the new folder
fdr_df.to_csv(os.path.join(results_folder, "GSEA_overlap_fdr_values.csv"), index=False)
power_df.to_csv(os.path.join(results_folder, "GSEA_overlap_power_values.csv"), index=False)

    