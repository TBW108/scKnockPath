
import os
import scanpy as sc
import numpy as np
import sys
sys.path.append("./model")
from model_lambda import scKnockPath
from groupyr import LogisticSGLCV
from sklearn.preprocessing import StandardScaler
import joblib


folder_path='./results/simulation_exp/ablation_sample'
    
seeds = np.arange(46, 76)
n_samples = np.linspace(600, 1800, 6, dtype=int)

for n_sample in n_samples:
    print(f'n_sample={n_sample}')
    fdrs=[]
    powers=[]
    for seed in seeds:
        adata=sc.read_h5ad(f'./data/simulation/sample_data/n_sample_beta{n_sample}/simu_scRNAseq_100pathways_{n_sample}sample_seed={seed}.h5ad')

        # # preprocess data
        # sc.pp.normalize_total(adata)
        # sc.pp.log1p(adata)
        adata.layers['lognorm']=adata.X.copy()
        true_effect_pathways=adata.uns['effect_pathways']
        geneset_dict=adata.uns['all_pathways']
        obs_y='cell_type'
        gene_thresh=1
        
        # create a folder to store results
        subfolder_name = f"n_sample={n_sample}" 
        subfolder_path = os.path.join(folder_path, subfolder_name)  
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)
            print("sub文件夹已创建！")
            
        print(f'data_seed={seed}')
        # fit the model
        model=scKnockPath()
        X, Xh, y, groups_list = model.prepare_data(
            adata=adata,
            layer="lognorm",
            obs_y=obs_y,
            genesets=geneset_dict,
            gene_thresh=gene_thresh,
            class1=1.0,
            class2=0.0,
        )
        Xh = StandardScaler().fit_transform(Xh)

        alphas=np.logspace(-1, -3, num=10)
        print(f'start fitting')
        sgl_model=LogisticSGLCV(l1_ratio=0,groups=groups_list,alphas=alphas,cv=3,n_jobs=-1)
        sgl_model.fit(Xh,y)
        
        final_pathways=model.pathway_names[sgl_model.chosen_groups_]
        
        joblib.dump(sgl_model, os.path.join(subfolder_path, f'sgl_model_seed_{seed}.pkl'))
        
        real_fdr=np.logical_not(np.isin(final_pathways,true_effect_pathways)).sum()/len(final_pathways)
        real_power=np.isin(true_effect_pathways,final_pathways).sum()/len(true_effect_pathways)
        print(f'OGL selects {len(final_pathways)}pathways')
        print(final_pathways,f'real_power={real_power},real_fdr={real_fdr}')
        fdrs.append(real_fdr)
        powers.append(real_power)

        print('-----------------')
        
    print(f'mean_fdrs={fdrs}')
    print(f'mean_powers={powers}')
    np.savez(f'{subfolder_path}/n_sample={n_sample}_results.npz',fdrs=fdrs,powers=powers)
    print('===============================================================================')