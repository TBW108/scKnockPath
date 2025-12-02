import os
import scanpy as sc
import numpy as np
import knockpy
import sys
sys.path.append("./model")
from util import plot_correlation  # 更简洁的导入方式
# from model_GMM import scKnockPath  # 更简洁的导入方式
from model_lambda import scKnockPath

knockoff_seed=42

# create a folder to store results
folder_name = f"{__file__.split('/')[-1].split('.')[0]}" 
folder_path = os.path.join('./results/simulation_exp', folder_name)  
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    print("文件夹已创建！")
    

n_overlap=4
seeds=np.arange(46,76)
n_sparsitys=np.arange(6,17,2)

for n_sparsity in n_sparsitys:
    print(f'n_sparsity={n_sparsity}')
    fdrs=[]
    powers=[]
    for seed in seeds:
        adata=sc.read_h5ad(f'./data/simulation/sparsity_data/n_sparsity_beta{n_sparsity}/simu_scRNAseq_100pathways_{n_sparsity}sparsity_seed={seed}.h5ad')

        # preprocess data
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.layers['lognorm']=adata.X.copy()
        true_effect_pathways=adata.uns['effect_pathways']
        geneset_dict=adata.uns['all_pathways']
        obs_y='cell_type'
        gene_thresh=1
        
        # create a folder to store results
        subfolder_name = f"n_sparsity={n_sparsity}" 
        subfolder_path = os.path.join(folder_path, subfolder_name)  
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)
            print("sub文件夹已创建！")
        
        # generate knockoff
        print(f'data_seed={seed}')
        np.random.seed(knockoff_seed)
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
        
        load=True
        datafile_name=f'100p_{n_sparsity}sparsity_data_seed={seed}.npz'
        data_file = os.path.join(subfolder_path, datafile_name)
        if os.path.exists(data_file) and load:
            data = np.load(data_file,allow_pickle=True)
            Xc = data['Xc']
            y = data['y']
        else:
            print('saving data')
            # Xc,X,Xk = model.knockoff_sampler(X,Xh,method='GMM')
            Xc,X,Xk = model.knockoff_sampler(X,Xh,method='sdp')
            np.savez(data_file, Xc=Xc, y=y)
        
        lambdas=np.logspace(-1, -3, num=15)
        
        'alpha one by one test'
        fdr = 0.2
        model.fit(Xc, y, groups_list, alphas=lambdas,thresh=0.01)
        final_selected_pathways = model.select_pathway(fdr=fdr, offset=0)
        
        model.save_model(f'{subfolder_path}/100pathways_{n_sparsity}sparsity_seed={seed}_lambda.pkl')
        
        selected_pathways = final_selected_pathways
        real_fdr = np.logical_not(
            np.isin(selected_pathways, true_effect_pathways)
        ).sum() / len(selected_pathways)
        real_power = np.isin(true_effect_pathways, selected_pathways).sum() / len(
            true_effect_pathways
        )
        print(
            f"target_fdr={fdr}, scKnockPath selected {len(selected_pathways)}pathways"
        )
        print(selected_pathways, f"real_power={real_power},real_fdr={real_fdr}")
        fdrs.append(real_fdr)
        powers.append(real_power)
        print("-" * 40)
        
        # for fdr in [0.2]:
        #     selected_pathways=model.select_pathway(fdr=fdr,offset=0)
        #     real_fdr=np.logical_not(np.isin(selected_pathways,true_effect_pathways)).sum()/len(selected_pathways)
        #     real_power=np.isin(true_effect_pathways,selected_pathways).sum()/len(true_effect_pathways)
        #     print(f'target_fdr={fdr}, scKnockPath selected {len(selected_pathways)}pathways')
        #     print(selected_pathways,f'real_power={real_power},real_fdr={real_fdr}')
        #     fdrs.append(real_fdr)
        #     powers.append(real_power)
        
        'SS_selection'
        # fdr=0.2
        # model.SS_fit(Xc,y,groups_list,alphas,fdr=0.2,N=10,N_sample=400)
        
        # model.save_model(f'{subfolder_path}/100pathways_{n_sparsity}sparsity_seed={seed}_SS.pkl')
        # final_pathways=model.SS_select_pathway(th=0.3)
        
        # real_fdr=np.logical_not(np.isin(final_pathways,true_effect_pathways)).sum()/len(final_pathways)
        # real_power=np.isin(true_effect_pathways,final_pathways).sum()/len(true_effect_pathways)
        # print(f'target_fdr={fdr}, scKnockPath selected {len(final_pathways)}pathways')
        # print(final_pathways,f'real_power={real_power},real_fdr={real_fdr}')
        # fdrs.append(real_fdr)
        # powers.append(real_power)
        
        'CV_selection'
        # model.fit(Xc,y,groups_list,n_jobs=-1,useCV=True,alphas=alphas,cv=4)
        # model.save_model(f'{subfolder_path}/100pathways_{n_overlap}overlapadd_seed={seed}_CV.pkl')
        # for fdr in [0.2]:
        #     final_pathways=model.select_pathway(fdr=fdr,offset=0)
        #     real_fdr=np.logical_not(np.isin(final_pathways,true_effect_pathways)).sum()/len(final_pathways)
        #     real_power=np.isin(true_effect_pathways,final_pathways).sum()/len(true_effect_pathways)
        #     print(f'fdr={fdr}, scKnockPath selected {len(final_pathways)}pathways')
        #     print(f'real_power={real_power},real_fdr={real_fdr}')
        #     fdrs.append(real_fdr)
        #     powers.append(real_power)
            
            
        print('-----------------')
        
    print(f'mean_fdrs={fdrs}')
    print(f'mean_powers={powers}')

    results_folder_path = os.path.join(folder_path, "results")
    
    if not os.path.exists(results_folder_path):
        os.makedirs(results_folder_path)
        print("results文件夹已创建！")
    np.savez(f'{results_folder_path}/n_sparsity={n_sparsity}_results.npz',fdrs=fdrs,powers=powers)
    print('===============================================================================')