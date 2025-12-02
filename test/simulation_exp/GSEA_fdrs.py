import numpy as np
import pandas as pd
import gseapy as gp
import scanpy as sc
import os
import warnings


def run_gsea(filename, obs_key="cell_type", seed=42):
    adata = sc.read_h5ad(filename)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    df = pd.DataFrame(adata.layers["lognorm"], columns=adata.var_names)
    geneset_dict = adata.uns["all_pathways"]
    effect_pathways = np.array(adata.uns["effect_pathways"], dtype=str)

    res = gp.gsea(
        data=df.T,
        gene_sets={
            pathway: genes.astype(str).tolist()
            for pathway, genes in geneset_dict.items()
        },
        cls=adata.obs[obs_key],
        permutation_num=1000,
        method="t_test",
        min_size=1,
        seed=seed,
    )
    return res.res2d, effect_pathways


def summarize_fdr(res_df, true_effect_pathways, fdr_thresholds):
    real_fdrs = []
    real_powers = []

    for fdr in fdr_thresholds:
        selected = res_df.loc[res_df["FDR q-val"] <= fdr, "Term"].to_numpy(dtype=str)
        selected_len = len(selected)

        if selected_len == 0:
            real_fdrs.append(0.0)
            real_powers.append(0.0)
            print(f"There are 0 pathways selected with FDR<={fdr}")
            continue

        false_discoveries = np.count_nonzero(
            ~np.isin(selected, true_effect_pathways)
        )
        true_discoveries = np.count_nonzero(
            np.isin(true_effect_pathways, selected)
        )

        real_fdrs.append(false_discoveries / selected_len)
        real_powers.append(true_discoveries / len(true_effect_pathways))
        print(f"There are {selected_len} pathways selected with FDR<={fdr}")

    return real_fdrs, real_powers


folder_name = f"{__file__.split('/')[-1].split('.')[0]}"
folder_path = "./results/simulation_exp/GSEA_fdrs_results"
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    print("文件夹已创建！")

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

fdr_values = []
power_values = []
fdr_thresholds = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3]
n_overlap = 4

path = f"./data/simulation/overlap_data/n_overlap_beta{n_overlap}"
for seed in range(46, 76):
    print(f"seed={seed}")
    res_df, effect_pathways = run_gsea(
        f"{path}/simu_scRNAseq_100pathways_{n_overlap}overlaplog_seed={seed}.h5ad"
    )
    fdr, power = summarize_fdr(res_df, effect_pathways, fdr_thresholds)
    fdr_values.append(fdr)
    power_values.append(power)

    print(f"real_fdr={fdr}, real_power={power}")
    print("--------------------------------------------------------------------")
    print(fdr_values)
    print(power_values)

    print("===============================================")

fdr_df = pd.DataFrame(fdr_values,columns=[str(fdr) for fdr in fdr_thresholds])
power_df = pd.DataFrame(power_values,columns=[str(fdr) for fdr in fdr_thresholds])


fdr_df.to_csv(os.path.join(folder_path, "GSEAfdrs_fdr_values.csv"), index=False)
power_df.to_csv(os.path.join(folder_path, "GSEAfdrs_power_values.csv"), index=False)
