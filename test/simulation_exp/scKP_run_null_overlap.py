import os
import scanpy as sc
import numpy as np
import knockpy
import seaborn as sns
import pandas as pd
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
    

n_overlaps = np.arange(0, 11, 2)
# seeds = np.arange(46, 56)
seeds = np.arange(46, 76)

for n_overlap in n_overlaps:
    print(f"n_overlap={n_overlap}")
    n_selected_pathways = []
    for seed in seeds:
        # load data
        adata = sc.read_h5ad(
            f"./data/simulation/null_overlap_data/n_overlap_beta{n_overlap}/simu_scRNAseq_100pathways_{n_overlap}overlap_null_seed={seed}.h5ad"
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

        # create a folder to store results
        subfolder_name = f"n_overlap={n_overlap}"
        subfolder_path = os.path.join(folder_path, subfolder_name)
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)
            print("sub文件夹已创建！")

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

        load = False
        datafile_name = f"100p_{n_overlap}overlap_knockoffdata_seed={seed}.npz"
        data_file = os.path.join(subfolder_path, datafile_name)
        if os.path.exists(data_file) and load:
            data = np.load(data_file, allow_pickle=True)
            Xc = data["Xc"]
            y = data["y"]
        else:
            print("saving data")
            Xc, X, Xk = model.knockoff_sampler(X, Xh, method="sdp")
            np.savez(data_file, Xc=Xc, y=y)

        # plot_correlation(X, Xk, subfolder_path, seed)
        print(f"shape of Xc={Xc.shape},shape of y={y.shape}")
        lambdas = np.logspace(-1, -3, num=15)

        "alpha one by one test"
        fdr = 0.2
        model.fit(Xc, y, groups_list, alphas=lambdas,thresh=-0.001)
        model.save_model(
            f"{subfolder_path}/100pathways_{n_overlap}overlap_null_seed={seed}_lambda.pkl"
        )
        final_selected_pathways = model.select_pathway(fdr=fdr, offset=0)
        n_selected_pathways.append(len(final_selected_pathways))
        print(f"n_selected_pathways=len(final_selected_pathways)={len(final_selected_pathways)}")
        print("-" * 40)


    # store the result
    print(f"mean_n_selected_pathways={np.mean(n_selected_pathways)}")

    # Create a new results folder under folder_path
    results_folder_path = os.path.join(folder_path, "results")
    if not os.path.exists(results_folder_path):
        os.makedirs(results_folder_path)
        print("results文件夹已创建！")

    # Save results in the new results folder
    results_file_path = os.path.join(
        results_folder_path, f"n_overlap={n_overlap}_results.npz"
    )
    np.savez(results_file_path, n_selected_pathways=n_selected_pathways)
    print(
        "==============================================================================="
    )
