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
    fdr = 0.2
    result = np.array(res.res2d[(res.res2d['FDR q-val'] <= fdr)]['Term'])
    
    real_fdr = np.logical_not(np.isin(result, true_effect_pathways)).sum() / len(result)
    real_power = np.isin(true_effect_pathways, result).sum() / len(true_effect_pathways)
    print(f'real_fdr={real_fdr}')
    print(f'real_power={real_power}')
    print(f'There are {len(result)} pathways selected with FDR<={fdr}')
    
    return real_fdr, real_power, len(result)



warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

fdr_list = {}
power_list = {}
n_selected_pathways_list = {}

for beta_value in [0, 0.5, 1, 1.5, 2, 3]:
    print(f'beta={beta_value}')
    path = f"./data/simulation/beta_data/beta={beta_value}"
    fdr_values = []
    power_values = []
    n_selected_pathways_values = []
    for seed in range(46, 76):
        print(f'seed={seed}')
        fdr, power, n_selected = GSEA_test(f"{path}/simu_scRNAseq_100pathways_beta={beta_value}_seed={seed}.h5ad")
        fdr_values.append(fdr)
        power_values.append(power)
        n_selected_pathways_values.append(n_selected)
        print(f'real_fdr={fdr}, real_power={power}, n_selected_pathways={n_selected}')
        print('--------------------------------------------------------------------')
            
    print(fdr_values)
    print(power_values)
    print(n_selected_pathways_values)
    
    fdr_list[str(beta_value)] = fdr_values
    power_list[str(beta_value)] = power_values
    n_selected_pathways_list[str(beta_value)] = n_selected_pathways_values
    print('===============================================')

fdr_df = pd.DataFrame(fdr_list)
power_df = pd.DataFrame(power_list)
n_selected_pathways_df = pd.DataFrame(n_selected_pathways_list)
# Create a folder to store results
results_folder = "./results/simulation_exp/GSEA_beta_results"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("GSEA_beta_results folder created!")

# Save results in the new folder
fdr_df.to_csv(os.path.join(results_folder, "GSEA_beta_fdr_values.csv"), index=False)
power_df.to_csv(os.path.join(results_folder, "GSEA_beta_power_values.csv"), index=False)
n_selected_pathways_df.to_csv(os.path.join(results_folder, "GSEA_beta_n_selected_values.csv"), index=False)

    