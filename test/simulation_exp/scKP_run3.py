import os
import scanpy as sc
import numpy as np
import sys
sys.path.append("./model")
from model_lambda import scKnockPath
np.random.seed(42)
# create a folder to store results
folder_name = f"{__file__.split('/')[-1].split('.')[0]}"
folder_path = os.path.join("./results/simulation_exp", folder_name)
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    print("文件夹已创建！")

seeds = np.arange(46, 76)
n_samples = np.linspace(600, 1800, 6, dtype=int)



for n_sample in n_samples:
    print(f"n_samples={n_samples}")
    fdrs = []
    powers = []
    for seed in seeds:
        # load data
        adata = sc.read_h5ad(
            f"./data/simulation/sample_data/n_sample_beta{n_sample}/simu_scRNAseq_100pathways_{n_sample}sample_seed={seed}.h5ad"
        )
        adata.layers["lognorm"] = adata.X.copy()
        true_effect_pathways = adata.uns["effect_pathways"]
        geneset_dict = adata.uns["all_pathways"]
        obs_y = "cell_type"
        gene_thresh = 1

        # create a folder to store results
        subfolder_name = f"n_sample={n_sample}"
        subfolder_path = os.path.join(folder_path, subfolder_name)
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)
            print("sub文件夹已创建！")

        # generate knockoff
        print(f"data_seed={seed}")
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
        datafile_name = f"100p_{n_sample}samples_data_seed={seed}.npz"
        data_file = os.path.join(subfolder_path, datafile_name)
        if os.path.exists(data_file) and load:
            data = np.load(data_file, allow_pickle=True)
            Xc = data["Xc"]
            y = data["y"]
        else:
            print("saving data")
            # Xc,X,Xk = model.knockoff_sampler(X,Xh,method='GMM')
            Xc, X, Xk = model.knockoff_sampler(X, Xh, method="sdp")
            np.savez(data_file, Xc=Xc, y=y)

        lambdas = np.logspace(-1, -3, num=15)
        fdr = 0.2
        model.fit(Xc, y, groups_list, alphas=lambdas)

        model.save_model(f'{subfolder_path}/100pathways_{n_sample}samples_seed={seed}_lambda.pkl')
        final_selected_pathways = model.select_pathway(fdr=fdr, offset=0)

        selected_pathways = final_selected_pathways
        # selected_pathways=model.select_pathway(fdr=fdr,offset=0)
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

        print("-----------------")

    print(f"mean_fdrs={fdrs}")
    print(f"mean_powers={powers}")

    results_folder_path = os.path.join(folder_path, "results")
    if not os.path.exists(results_folder_path):
        os.makedirs(results_folder_path)
        print("results文件夹已创建！")
    np.savez(
        f"{results_folder_path}/n_samples={n_sample}_results.npz",
        fdrs=fdrs,
        powers=powers,
    )
    print(
        "==============================================================================="
    )
