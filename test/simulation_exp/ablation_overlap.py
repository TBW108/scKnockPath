
import os
import scanpy as sc
import numpy as np
import sys
sys.path.append("./model")
from model_lambda import scKnockPath
from groupyr import LogisticSGLCV,LogisticSGL
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import joblib
np.random.seed(42)

n_overlaps=np.arange(0,11,2)
seeds=np.arange(46,76)
use_CV=True
cv=3 if use_CV else None
thresh=0.01

folder_path='./results/simulation_exp/ablation_overlap'
for n_overlap in n_overlaps:
    print(f'n_overlap={n_overlap}')
    fdrs=[]
    powers=[]
    for seed in seeds:
        adata=sc.read_h5ad(f'./data/simulation/overlap_data/n_overlap_beta{n_overlap}/simu_scRNAseq_100pathways_{n_overlap}overlaplog_seed={seed}.h5ad')

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
        subfolder_name = f"n_overlap={n_overlap}" 
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


        alphas=np.logspace(-1, -3, num=10)
        print(f'start fitting')
        Xh = StandardScaler().fit_transform(Xh)
        
        best_model=LogisticSGLCV(l1_ratio=0,groups=groups_list,alphas=alphas,cv=cv,n_jobs=-1)
        best_model.fit(Xh,y)
        print(f'fit done. It selects {len(best_model.chosen_groups_)} pathways')
        
        # if use_CV:
        #     for i,alpha in enumerate(alphas):
        #         print(f'alpha={alpha}')
        #         # Use KFold cross-validation to fit Xc and y with the given alpha
        #         kf = KFold(n_splits=cv if isinstance(cv, int) else 5, shuffle=True, random_state=42)
        #         scores = []
        #         sgl_model = LogisticSGL(l1_ratio=0,groups=groups_list,alpha=alpha)
        #         for train_index, test_index in kf.split(Xh):
        #             X_train, X_test = Xh[train_index], Xh[test_index]
        #             y_train, y_test = y[train_index], y[test_index]
        #             sgl_model.fit(X_train, y_train)
        #             scores.append(sgl_model.score(X_test, y_test))
        #         print(f"KFold CV mean score: {np.mean(scores):.4f}")
        #         score=np.mean(scores)
                
        #         # test if there are pathways are selected   
        #         sgl_model.fit(Xh,y)
        #         selected_pathways=model.pathway_names[sgl_model.chosen_groups_]
        #         print(f'OGL selects {len(selected_pathways)} pathways')
        #         print(selected_pathways)
        #         # when the score does not increase, break
        #         if (score-best_score)<thresh and len(selected_pathways) != 0:  
        #         # Threshold for early dropping
        #             print(f'Early drop at alpha={alpha}, score={score} drops or keeps constant.')
        #             break
                    
        #         # update the best score and best model
        #         best_model=sgl_model if best_model is None or score > best_score else best_model
        #         best_score=score if score >best_score else best_score
        #         print("*"*40)
        
        joblib.dump(best_model, os.path.join(subfolder_path, f'sgl_model_seed_{seed}.pkl'))
        final_pathways=model.pathway_names[best_model.chosen_groups_]
        real_fdr=np.logical_not(np.isin(final_pathways,true_effect_pathways)).sum()/len(final_pathways)
        real_power=np.isin(true_effect_pathways,final_pathways).sum()/len(true_effect_pathways)
        print(final_pathways,f'real_power={real_power},real_fdr={real_fdr}')
        fdrs.append(real_fdr)
        powers.append(real_power)
        

        print('-----------------')
        
    print(f'mean_fdrs={fdrs}')
    print(f'mean_powers={powers}')

    np.savez(f'{subfolder_path}/n_overlap={n_overlap}_results.npz',fdrs=fdrs,powers=powers)
    print('===============================================================================')