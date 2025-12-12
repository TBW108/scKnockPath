import os
import scanpy as sc
import numpy as np
import sys
import pickle as pkl
# 添加上级目录
sys.path.append("./model")
from model_lambda import scKnockPath

def get_fdr_power(selected_pathways, true_effect_pathways, beta_value):
    if len(selected_pathways) == 0:
        real_fdr = 0.0
        if beta_value == 0:
            real_power = None # beta=0 时应该记录 0 的次数
        else:
            real_power = 0
    else:
        if beta_value == 0:
            real_fdr = 1.0
            real_power = None
        else:
            real_fdr = np.logical_not(
                np.isin(selected_pathways, true_effect_pathways)
            ).sum() / len(selected_pathways)
            real_power = np.isin(true_effect_pathways, selected_pathways).sum() / len(true_effect_pathways)
    return real_fdr, real_power

knockoff_seed = 42

# create a folder to store results
folder_name = f"{__file__.split('/')[-1].split('.')[0]}"
folder_path = os.path.join("./results/simulation_exp", folder_name)
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    print("文件夹已创建！")
    

beta_values = [0, 0.5, 1, 1.5, 2, 3]
n_overlap = 4
seeds = np.arange(46, 76)
fdr = 0.2
for beta_value in beta_values:
    print(f"beta_value={beta_value}")
    fdrs = []
    powers = []
    
    # create a folder to store results
    subfolder_name = f"beta={beta_value}"
    subfolder_path = os.path.join(folder_path, subfolder_name)
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)
        print("sub文件夹已创建！")
        
    for seed in seeds:
        # load data
        adata = sc.read_h5ad(
            f"./data/simulation/beta_data/beta={beta_value}/simu_scRNAseq_100pathways_beta={beta_value}_seed={seed}.h5ad"
        )       
        
        #preprocess data
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=10)
        
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.layers["lognorm"] = adata.X
        true_effect_pathways = adata.uns["effect_pathways"]
        geneset_dict = adata.uns["all_pathways"]
        obs_y = "cell_type"
        gene_thresh = 1


        print(f"data_seed={seed}")
        
        # if the pkl file already exists, skip
        pkl_file_path = f"./results/simulation_exp/scKP_run_beta/beta={beta_value}/100pathways_beta={beta_value}_seed={seed}_lambda.pkl"
        if os.path.exists(pkl_file_path):
            print("load model from pkl file")
            model=pkl.load(open(pkl_file_path,'rb'))
            selected_pathways = model.select_pathway(fdr=fdr, offset=0)
            
            # store the result
            real_fdr, real_power = get_fdr_power(
            selected_pathways, true_effect_pathways, beta_value
            )
            print(
                f"target_fdr={fdr}, scKnockPath selected {len(selected_pathways)}pathways" )
            print(f"real_power={real_power},real_fdr={real_fdr}")
            fdrs.append(real_fdr)
            powers.append(real_power)
            print("-" * 40)
            continue
        
        # establish the model
        print("establish the model")
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
        datafile_name = f"100p_{n_overlap}overlap_knockoffdata_seed={seed}.npz"
        temp_path = f'./results/simulation_exp/scKP_run/n_overlap={n_overlap}'
        data_file = os.path.join(temp_path, datafile_name)
        
        data_loaded = False
        if os.path.exists(data_file) and load:
            print("loading data")
            try:
                data = np.load(data_file, allow_pickle=True)
                Xc = data["Xc"]
                y = np.array(adata.obs[obs_y])
                data_loaded = True
            except Exception as e:
                print(f"Failed to load data from {data_file}: {e}")
                print("Regenerating data...")

        if not data_loaded:
            print("saving data")
            np.random.seed(knockoff_seed)
            Xc, X, Xk = model.knockoff_sampler(X, Xh, method="sdp")
            np.savez(data_file, Xc=Xc, y=y)

        print(f"shape of Xc={Xc.shape},shape of y={y.shape}")
        lambdas = np.logspace(-1, -3, num=15)

        "alpha one by one test"

        model.fit(Xc, y, groups_list, alphas=lambdas,thresh=0.01)
        model.save_model(
            f"{subfolder_path}/100pathways_beta={beta_value}_seed={seed}_lambda.pkl"
        )
        final_selected_pathways = model.select_pathway(fdr=fdr, offset=0)

        selected_pathways = final_selected_pathways
        # 如果选择的通路为空，beta值也为 0，则 true_effect_pathways 也为空，此时定义 fdr=0 和 power 都为 1. 如果 beta 为非0值，则 power 定义为0，fdr定义为0
        real_fdr, real_power = get_fdr_power(
            selected_pathways, true_effect_pathways, beta_value
        )
        print(
            f"target_fdr={fdr}, scKnockPath selected {len(selected_pathways)}pathways"
        )
        print(f"real_power={real_power},real_fdr={real_fdr}")
        fdrs.append(real_fdr)
        powers.append(real_power)
        print("-" * 40)

    # store the result
    print(f"mean_fdrs={fdrs}")
    print(f"mean_powers={powers}")

    # Create a new results folder under folder_path
    results_folder_path = os.path.join(folder_path, "results")
    if not os.path.exists(results_folder_path):
        os.makedirs(results_folder_path)
        print("results文件夹已创建！")

    # Save results in the new results folder
    results_file_path = os.path.join(
        results_folder_path, f"beta={beta_value}_results.npz"
    )
    np.savez(results_file_path, fdrs=fdrs, powers=powers)
    print(
        "==============================================================================="
    )
