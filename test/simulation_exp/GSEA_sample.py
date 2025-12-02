'''
Based on the rank of pathways using JI, test whether scKnockPath, OGLasso with CV, SCPA, GSEA can select true pathways
'''
import numpy as np
import pandas as pd
import gseapy as gp
from groupyr import LogisticSGL,LogisticSGLCV
import scanpy as sc
import pandas as pd
import pickle
import os
# --------------------------------------------------------------------------------

def GSEA_test(filename):
    adata = sc.read_h5ad(filename)
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()

    true_effect_pathways = adata.uns['effect_pathways']

    df = pd.DataFrame(adata.layers['lognorm'], columns=adata.var_names)
    geneset_dict = adata.uns['all_pathways']
    obs_y = 'cell_type'

    res = gp.gsea(
        data=df.T,
        gene_sets={pathway: genes.astype(str).tolist() for pathway, genes in geneset_dict.items()},
        cls=adata.obs[obs_y],
        permutation_num=2000,
        method='t_test',
        min_size=1,
        seed=42
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


fdr_list = {}
power_list = {}

for n_sample in np.linspace(600,1800,6):
    n_sample = int(n_sample)
    print(f'n_sample={n_sample}')
    path = f"./data/simulation/sample_data/n_sample_beta{n_sample}"
    fdr_values = []
    power_values = []
    
    for seed in range(46, 76):
        print(f'seed={seed}')
        fdr, power = GSEA_test(f"{path}/simu_scRNAseq_100pathways_{n_sample}sample_seed={seed}.h5ad")
        fdr_values.append(fdr)
        power_values.append(power)
        print(f'real_fdr={fdr}, real_power={power}')
        print('--------------------------------------------------------------------')
            
    print(fdr_values)
    print(power_values)
    
    fdr_list[str(n_sample)] = fdr_values
    power_list[str(n_sample)] = power_values
    print('===============================================')

fdr_df = pd.DataFrame(fdr_list)
power_df = pd.DataFrame(power_list)
# Create a folder to store results
results_folder = "./results/simulation_exp/GSEA_sample_results"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("GSEA_sample_results folder created!")

# Save results in the new folder
fdr_df.to_csv(os.path.join(results_folder, "GSEA_sample_fdr_values.csv"), index=False)
power_df.to_csv(os.path.join(results_folder, "GSEA_sample_power_values.csv"), index=False)

fdr_df.to_csv("./results/simulation_exp/GSEA_fdr_values.csv", index=False)
power_df.to_csv("./results/simulation_exp/GSEA_power_values.csv", index=False)


# for i in range(0, 11):
#     print(f'n_overlap={i}')
#     GSEA_test(f'/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_100pathways_{i}overlap.h5ad')

# n_sparsity = np.arange(10, 21, 2)
# for i in n_sparsity:
#     print(f'n_sparse={i}')
#     GSEA_test(f'/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_100pathways_4overlap{i}.h5ad')

# adata=sc.read_h5ad('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_5pathways_overlapadd.h5ad')

# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata)
# adata.layers['lognorm']=adata.X.copy()

# true_effect_pathways=adata.uns['effect_pathways']

# df=pd.DataFrame(adata.layers['lognorm'],columns=adata.var_names)
# geneset_dict=adata.uns['all_pathways']
# obs_y='cell_type'
# # ground_truth=adata.uns['ground_truth']

# res = gp.gsea(
#     data=df.T,
#     gene_sets={pathway:genes.astype(str).tolist() for pathway,genes in geneset_dict.items()},  # 直接传递字典
#     cls=adata.obs[obs_y],
#     permutation_num=1000,  # 设置置换次数
#     method='t_test',  # 排序方法，可以选择'signal_to_noise'等
#     min_size=1# 5/15is another option
#     )

# fdrs=[0.2,0.1,0.05]
# print(res.res2d)
# for fdr in fdrs:
#     result=np.array(res.res2d[(res.res2d['FDR q-val'] <= fdr)]['Term'])
    
#     print(f'There are {len(result)} pathways selected with FDR<={fdr}')
#     print(np.isin(true_effect_pathways,result))
    # print(result)


# sig_pathways.append(np.array(res.res2d[(res.res2d['FDR q-val'] <= 0.2)]['Term']))
# sig_pathways.append(np.array(res.res2d[(res.res2d['FDR q-val'] <= 0.1)]['Term']))
# sig_pathways.append(np.array(res.res2d[(res.res2d['FDR q-val'] <= 0.05)]['Term']))



# selected_id=np.zeros(len(geneset_dict),dtype=bool)
# genesets_names=np.array(list(geneset_dict.keys()))
# selected_id[np.isin(genesets_names,sig_pathways)]=True        
# real_fdr,real_power=get_fdr_power(selected_id,ground_truth)

# print(f"{model_name}:{len(sig_pathways)} pathways are selected, real_fdr={real_fdr}, real power={real_power}")






# # 计算每个不同JI的pathway能否被选出来
# for method in methods:
#     # create the folder of the method to store results
#     subfolder_path = os.path.join(folder_path, method)  # the path that stores results for the method
#     if not os.path.exists(subfolder_path):
#         os.makedirs(subfolder_path)  
#         print('Create the folder success!')
    
#     for pathway in pathway_JI['pathways']:
#         # read in the adata
#         adata = sc.read_h5ad(f'scKnockPath/compare_results/simulation/pathwayJI/simu_scRNAseq_{pathway}.h5ad')
#         sc.pp.normalize_total(adata)
#         sc.pp.log1p(adata)
#         adata.layers['lognorm']=adata.X.copy()
        
#         JI=adata.uns['JI']
#         true_effect_pathways=adata.uns['effect_pathways']
#         print(adata.uns['JI'])
#         print(adata.uns['effect_pathways'])
#         real_fdr,real_power,selected_pathways=compare(adata,'cell_type',
#                                                     adata.uns['pathways'],
#                                                     method,
#                                                     adata.uns['ground_truth']) 
        
#         sim_res=dict()
#         sim_res['real_fdr']=real_fdr
#         sim_res['real_power']=real_power
#         sim_res['selected_pathways']=selected_pathways
#         sim_res['true_effect_pathways']=true_effect_pathways
        
#         with open(f'{subfolder_path}/path={pathway}.pkl', 'wb') as f:
#             pickle.dump(sim_res, f)
        

    

