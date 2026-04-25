import os
import sys
import numpy as np
import scanpy as sc
import warnings
import time
sys.path.append("./model")
# from model_lambda_alt import scKnockPath
from sklearn.preprocessing import QuantileTransformer
from model_lambda import scKnockPath
from scipy.sparse import issparse
import argparse

# 过滤掉包含 'resource_tracker' 的警告
warnings.filterwarnings("ignore", message="resource_tracker: There appear to be .* leaked folder objects")

def rank_gauss_transform(adata, layer=None):
    """
    对单细胞矩阵的每个基因进行逆正态变换 (RankGauss)
    """
    # 提取表达矩阵 (假设是 log-norm 后的稠密或稀疏矩阵)
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X
        
    # 如果是稀疏矩阵，转换为稠密矩阵 (RankGauss 通常需要稠密计算)
    if hasattr(X, "toarray"):
        X = X.toarray()
        
    # 初始化 QuantileTransformer
    # output_distribution='normal' 是核心参数
    # n_quantiles 建议设为细胞数 (如果细胞数>1000，可设为1000以加速)
    n_cells = X.shape[0]
    qt = QuantileTransformer(
        n_quantiles=1000, 
        output_distribution='normal',
        random_state=42
    )
    
    # 对矩阵进行变换 (注意：fit_transform 是按列/按基因处理的，这正是我们需要的)
    X_gaussian = qt.fit_transform(X)
    
    # 将变换后的结果存回 adata
    adata.layers['rankgauss'] = X_gaussian
    return adata

parser = argparse.ArgumentParser(description="Run scKnockPath simulation experiments.")
parser.add_argument(
    "--sim_param",
    type=str,
    default="overlap",
    choices=["overlap", "sparsity", "sample", "signal",'nosignal'],
    help="Simulation experiment parameter to vary (default: overlap)",)

args = parser.parse_args()
param=args.sim_param

# 添加上级目录以导入自定义模型
if param =='overlap':
    pram_range=np.arange(0,11,5)
elif param =='sample':
    # pram_range=np.arange(600,1801,600)
    pram_range=[3000,4000,5000]
elif param =='signal':
    pram_range=[0.4,0.6,0.8]
elif param =='nosignal':
    pram_range=[0]

# 创建结果主文件夹
current_script_name = 'scKnockPath_' + param + '_results'
folder_path = os.path.join("./results/simulation_exp", current_script_name)
os.makedirs(folder_path, exist_ok=True)
print(f"主结果文件夹: {folder_path}")

# 结果子文件夹
results_folder_path = os.path.join(folder_path, "results")
os.makedirs(results_folder_path, exist_ok=True)

# 实验参数
data_seeds = np.arange(46, 56)
# data_seeds = [46]  # 或者 np.arange(46, 56) 如果需要多个
knockoff_seeds = [42]  # 或者 np.arange(42, 52) 如果需要多个
l1_ratio=0.3
knockoff_method="sdp"  # "sdp" 或 "permutation"
save_model=False # 保存模型
load_knockoff=True # 是否加载已生成的 knockoff 数据以节省时间
load_model=False
copula_rankgauss=True

for param_value in pram_range:
    print(f"\n{'='*30}\nProcessing n_{param}={param_value}\n{'='*30}")
    
    # 为当前 n_overlap 创建数据存储子文件夹
    if param=='nosignal':
        subfolder_name = f"n_signal={param_value}"
    else:
        subfolder_name = f"{param}={param_value}"
    subfolder_path = os.path.join(folder_path, subfolder_name)
    os.makedirs(subfolder_path, exist_ok=True)

    # 存储结果的列表
    fdrs = []
    powers = []
    n_selected_list = []
    # 内层循环：遍历数据集种子
    for data_seed in data_seeds:
        print(f"Processing data_seed={data_seed}...")
        
        # 1. 加载数据
        if param=='nosignal':
            data_path = f"./data/simulation/signal_data/signal={param_value}/simu_scRNAseq_100pathways_{param_value}signal_seed={data_seed}.h5ad"
        else:
            data_path = f"./data/simulation/{param}_data/{param}={param_value}/simu_scRNAseq_100pathways_{param_value}{param}_seed={data_seed}.h5ad"
        
        if not os.path.exists(data_path):
            print(f"Warning: Data file not found: {data_path}")
            # 如果文件确实，依然存入占位符（如 -1）以保证数组形状对齐，或者直接跳过（这会影响最终数组转换）
            # 这里选择打印警告并跳过当前种子的后续处理，但填充 NaN
            fdrs.append([np.nan] * len(knockoff_seeds))
            powers.append([np.nan] * len(knockoff_seeds))
            continue

        adata = sc.read_h5ad(data_path)
        # 2. 预处理
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        
        if copula_rankgauss:
            print("Applying RankGauss transformation...")
            adata=rank_gauss_transform(adata, layer=None) 
        else:
            print("Skipping RankGauss transformation, using log-normed data...")
            adata.layers['rankgauss'] = adata.X  # 如果不使用 copula rankgauss，直接复制原始 X
        print('processed data:')
        # print(adata.layers['rankgauss'][:10, :10])  # 打印前10行10列检查数据
        print('original data:')
        # print(adata.X[:10, :10])  # 打印原始数据前10行10列检查数据

        true_effect_pathways = adata.uns.get("effect_pathways", [])
        print(f"True effect pathways: {true_effect_pathways}")
        geneset_dict = adata.uns.get("all_pathways", {})
        obs_y = "cell_type"
        gene_thresh = 1

        repeat_fdrs = []
        repeat_powers = []
        repeat_n_selected = []

        # 3. 运行 Knockoff 分析
        for k_seed in knockoff_seeds:

            # 計算時間
            np.random.seed(k_seed)
            # 初始化模型并准备数据
            # 如果已经存在模型文件，可以选择加载以节省时间
            if copula_rankgauss:
                model_file_path = os.path.join(subfolder_path, 
                                           f"100pathways_{param_value}{param}_dataseed={data_seed}_knockoffseed={k_seed}_copula.pkl")
            else:
                model_file_path = os.path.join(subfolder_path, 
                                           f"100pathways_{param_value}{param}_dataseed={data_seed}_knockoffseed={k_seed}.pkl")
            
            if os.path.exists(model_file_path) and load_model:  
                print(f"Loading existing model from {model_file_path}")
                model = scKnockPath()
                model.load_model(model_file_path)
            else:
                model = scKnockPath()
                adata.obs[obs_y] = adata.obs[obs_y].astype(str)

                class1_val = "1"  # 对应原来的 1.0 或 1
                class2_val = "0"  # 对应原来的 0.0 或 0

                X, Xh, y, groups_list = model.prepare_data(
                    adata=adata,
                    layer="rankgauss",  # 使用 RankGauss 变换后的数据
                    obs_y=obs_y,
                    genesets=geneset_dict,
                    class1=class1_val, # 使用字符串 "1"
                    class2=class2_val, # 使用字符串 "0"
                )
                # l1_ratio=compute_data_driven_alpha(Xh,groups_list)
                start_time = time.time()
                
                # 加载或生成 Knockoff 数据
                if copula_rankgauss:
                    knockoff_file_name = f"100p_{param_value}{param}_dataseed={data_seed}_knockoffseed={k_seed}_copula.npz"
                else:
                    knockoff_file_name = f"100p_{param_value}{param}_dataseed={data_seed}_knockoffseed={k_seed}.npz"
                knockoff_file_path = os.path.join(subfolder_path, knockoff_file_name)
                
                if os.path.exists(knockoff_file_path) and load_knockoff:
                    print(f"Loading knockoff data from {knockoff_file_path}")
                    try:
                        data = np.load(knockoff_file_path, allow_pickle=True)
                        Xc = data["Xc"]
                        y_loaded = data["y"]
                        # 确保加载的 y 和 prepare_data 的 y 一致，或者直接使用加载的 y
                        if y_loaded.shape != y.shape:
                            print(f"Warning: Loaded y shape {y_loaded.shape} mismatch with current {y.shape}. Using generated Xc with current y.")
                            # 这是一个潜在风险点，如果样本没变，通常维度是一样的
                    except Exception as e:
                        print(f"Error loading knockoff data: {e}. Regenerating...")
                        Xc, X, Xk = model.knockoff_sampler(X, Xh, method=knockoff_method)
                        np.savez(knockoff_file_path, Xc=Xc, y=y)
                else:
                    print("Generating and saving knockoff data...")
                    # Xc, X, Xk = model.knockoff_sampler(X, Xh, method="sdp")
                    Xc, X, Xk = model.knockoff_sampler(X, Xh, method=knockoff_method)  
                    if save_model:
                        np.savez(knockoff_file_path, Xc=Xc, y=y)
                
                model.fit(Xc, y, groups_list,n_alphas=20,cv=3,l1_ratio=l1_ratio)
                y = y.astype(float)  # 确保 y 是数值类型以进行计算
                print(np.unique(y, return_counts=True))  # 检查 y 的类别分布
                # model.fit_W(Xc, y, groups_list)  # 直接计算 W 统计量，无需交叉验证拟合模型
                # model.fit(Xc, y, groups_list,cv=None,l1_ratio=0.1,alphas=[0.005])
                end_time = time.time()
                print(f"Model fitting took {end_time - start_time:.2f} seconds.")
                
                # 保存模型文件
                if copula_rankgauss:
                    model_save_path = os.path.join(subfolder_path, f"100pathways_{param_value}{param}_dataseed={data_seed}_knockoffseed={k_seed}_copula.pkl")
                else:
                    model_save_path = os.path.join(subfolder_path, f"100pathways_{param_value}{param}_dataseed={data_seed}_knockoffseed={k_seed}.pkl")
                
                if save_model:
                    model.save_model(model_save_path)
            
            # 选择通路并计算 FDR/Power
            target_fdr = 0.2
            final_selected_pathways = model.select_pathway(fdr=target_fdr, offset=0)
            
            def fdr_power(selected, true_effect):
                n_selected = len(selected)
                n_true = len(true_effect)
                
                if n_selected > 0:
                    # False Discoveries
                    false_discoveries = np.logical_not(np.isin(selected, true_effect)).sum()
                    real_fdr = false_discoveries / n_selected
                else:
                    real_fdr = 0.0 # 或者 np.nan，视你的统计需求而定，通常选不到即没有错误发现，FDR为0
                
                if n_true > 0:
                    true_positives = np.isin(true_effect, selected).sum()
                    real_power = true_positives / n_true
                else:
                    real_power = 0.0 # 理论上 true_effect_pathways 不应为空
                
                return real_fdr, real_power
            
            real_fdr, real_power = fdr_power(final_selected_pathways, true_effect_pathways)
            n_selected = len(final_selected_pathways)
                
            print(f"[KSeed={k_seed}] FDR={real_fdr:.4f}, Power={real_power:.4f}, Selected={n_selected}")
            
            
            
            final_selected_pathways_plus = model.select_pathway(fdr=target_fdr, offset=1)
            real_fdr_plus, real_power_plus = fdr_power(final_selected_pathways_plus, true_effect_pathways)
            n_selected = len(final_selected_pathways_plus)
            print(f"[KSeed={k_seed}] FDR+={real_fdr_plus:.4f}, Power={real_power_plus:.4f}, Selected={n_selected}")
            
            repeat_fdrs.append(real_fdr)
            repeat_powers.append(real_power)
            repeat_n_selected.append(n_selected)

        # 添加到总列表
        fdrs.append(repeat_fdrs)
        powers.append(repeat_powers)
        n_selected_list.append(repeat_n_selected)
        print("-" * 40)

    # 转换为 NumPy 数组以便保存
    fdrs_array = np.array(fdrs, dtype=float)
    powers_array = np.array(powers, dtype=float)
    n_selected_array = np.array(n_selected_list, dtype=int)
    print(f"all_FDRs {fdrs_array}, all_Powers {powers_array}, all_n_selected {n_selected_array}")
    
    if copula_rankgauss:
        print(f"Results for n_{param}={param_value} with Copula RankGauss:")
        if param=='nosignal':
            results_file_path = os.path.join(results_folder_path, f"n_signal={param_value}_results_copula.npz")
        else:
            results_file_path = os.path.join(results_folder_path, f"n_{param}={param_value}_results_copula.npz")
    else:
        print(f"Results for n_{param}={param_value} without Copula RankGauss:")
        if param=='nosignal':
            results_file_path = os.path.join(results_folder_path, f"n_signal={param_value}_results.npz")
        else:
            results_file_path = os.path.join(results_folder_path, f"n_{param}={param_value}_results.npz")
    
    # 保存为纯数值数组，无需 allow_pickle=True
    if save_model:
        np.savez(results_file_path, fdrs=fdrs_array, powers=powers_array, n_selected=n_selected_array)
        print(f"Saved results to {results_file_path}")

print("\nProcess Completed.")
