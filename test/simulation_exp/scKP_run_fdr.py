import os
import scanpy as sc
import numpy as np
import knockpy
import seaborn as sns
import pickle as pkl
import sys
# 添加上级目录
sys.path.append("./model")
# from util import plot_correlation  # 更简洁的导入方式
from model_lambda import scKnockPath
import matplotlib.pyplot as plt

knockoff_seed = 42

# create a folder to store results
folder_name = f"{__file__.split('/')[-1].split('.')[0]}"
folder_path = os.path.join("./results/simulation_exp", folder_name)
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    print("文件夹已创建！")
    

target_fdrs = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3]
seeds = np.arange(46, 47)  # 30 datasets
n_overlap = 4

# save the models
for seed in seeds:
    # load data
    adata = sc.read_h5ad(
        f"./data/simulation/overlap_data/n_overlap_beta{n_overlap}/simu_scRNAseq_100pathways_{n_overlap}overlaplog_seed={seed}.h5ad"
    )

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    # preprocess data
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    adata.layers["lognorm"] = adata.X.copy()
    true_effect_pathways = adata.uns["effect_pathways"]
    geneset_dict = adata.uns["all_pathways"]
    obs_y = "cell_type"
    gene_thresh = 1

    print(f"data_seed={seed}")
    np.random.seed(knockoff_seed)
    # fit the model
    model = scKnockPath()
    X, Xh, y, groups_list = model.prepare_data(
        adata=adata,
        layer="lognorm",
        obs_y=obs_y,
        genesets=geneset_dict,
        gene_thresh=gene_thresh,
        class1=1.0,
        class2=0.0,
    )

    load = True
    data_file = f"./results/simulation_exp/scKP_run/n_overlap={n_overlap}/100p_{n_overlap}overlap_knockoffdata_seed={seed}.npz"
    
    if os.path.exists(data_file) and load:
        data = np.load(data_file, allow_pickle=True)
        Xc = data["Xc"]
        y = data["y"]

    # plot_correlation(X, Xk, subfolder_path, seed)
    print(f"shape of Xc={Xc.shape},shape of y={y.shape}")
    lambdas = np.logspace(-1, -3, num=15)
    "alpha one by one test"
    model.fit(Xc, y, groups_list, alphas=lambdas,thresh=0.001)
    
    model.save_model(
        f"./results/simulation_exp/scKP_run_fdr/saved_models/overlap{n_overlap}_scKPmodel_seed={seed}.npz")
    



for fdr in target_fdrs:
    # record the fdr and power for each seed
    actual_fdrs = []
    actual_powers = []
    for seed in seeds:
        adata = sc.read_h5ad(
        f"./data/simulation/overlap_data/n_overlap_beta{n_overlap}/simu_scRNAseq_100pathways_{n_overlap}overlaplog_seed={seed}.h5ad"
            )
        true_effect_pathways = adata.uns["effect_pathways"]
        model=pkl.load(open(f"./results/simulation_exp/scKP_run_fdr/saved_models/overlap{n_overlap}_scKPmodel_seed={seed}.npz",'rb'))
        print(f"data_seed={seed}")

        final_selected_pathways = model.select_pathway(fdr=fdr, offset=0)
        selected_pathways = final_selected_pathways
        
        real_fdr = np.logical_not(
            np.isin(selected_pathways, true_effect_pathways)
        ).sum() / len(selected_pathways) if len(selected_pathways)>0 else 0
        
        real_power = np.isin(true_effect_pathways, selected_pathways).sum() / len(
            true_effect_pathways
        )
        print(
            f"target_fdr={fdr}, scKnockPath selected {len(selected_pathways)}pathways"
        )
        print(f"real_power={real_power},real_fdr={real_fdr}")
        actual_fdrs.append(real_fdr)
        actual_powers.append(real_power)
        print("-" * 40)
        
    # store the result
    print(f"mean_fdrs={actual_fdrs}")
    print(f"mean_powers={actual_powers}")

    # Save results in the new 
    results_file_path = os.path.join(
        folder_path,'results', f"fdr={fdr}_results.npz"
    )
    np.savez(results_file_path, fdrs=actual_fdrs, powers=actual_powers)
    print(
        "==============================================================================="
    )
