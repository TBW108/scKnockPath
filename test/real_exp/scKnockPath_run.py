import sys
import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
sys.path.append('./model/') # 添加model文件夹到系统路径
from model_lambda import scKnockPath
from sklearn.preprocessing import QuantileTransformer
import argparse
np.random.seed(42)

parser = argparse.ArgumentParser(description='scKnockPath pathway analysis')
parser.add_argument('--adata', type=str, required=True, help='Path to adata h5ad file')
parser.add_argument('--layer', type=str, default=None, help='Layer to use from adata')
parser.add_argument('--gene_names', type=str, default=None, help='Gene names obs (default: adata.var_names)')
parser.add_argument('--geneset_path', type=str, required=True, help='Path to pathway database file')
parser.add_argument('--l1_ratio', type=str, default='0.02', help='L1 ratio for elastic net regularization (e.g., 0.05)')
parser.add_argument('--obs_y', type=str, required=True, help='Column name for y variable in adata.obs')
parser.add_argument('--class1', type=str, required=True, help='First class name for y variable')
parser.add_argument('--class2', type=str, required=True, help='Second class name for y variable')
parser.add_argument('--fdr', type=float, default=0.2, help='FDR level for pathway selection (default: 0.2)')
parser.add_argument('--gene_thresh', type=int, default=15, help='Number of genes in each pathway')
parser.add_argument('--load_knockoff', type=bool, default=False, help='Whether to load knockoff data (default: False)')
parser.add_argument('--save_knockoff', type=bool, default=True, help='Whether to save knockoff data (default: True)')
parser.add_argument('--save_model', type=bool, default=True, help='Whether to save the trained model (default: True)')
parser.add_argument('--cell_type', type=str,default=None, help='Cell type for analysis')
parser.add_argument('--knockoff_seed', type=int, default=42, help='Random seed for knockoff generation (default: 42)')
parser.add_argument('--rank_gauss', type=bool, default=True, help='Whether to rank normalize the data (default: True)')
args = parser.parse_args()


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



# Extract filename prefix from adata path
adata_filename = os.path.basename(args.adata)
filename_prefix = os.path.splitext(adata_filename)[0]

# Create corresponding folder in ./data
data_folder = os.path.join('./results/real_exp', filename_prefix,'scKnockPath_results')
os.makedirs(data_folder, exist_ok=True)


adata=sc.read_h5ad(args.adata)

if args.cell_type is not None:
    if 'cell_type' not in adata.obs.columns:
        print("Error: 'cell_type' column not found in adata.obs. Please ensure your data contains cell type annotations.")
        sys.exit(1)
    adata = adata[adata.obs['cell_type'] == args.cell_type, :]
    print(f"Filtered data to {args.cell_type} cells: {adata.n_obs} cells remaining")
    
if args.rank_gauss:
    print("Applying RankGauss transformation to the data...")
    adata = rank_gauss_transform(adata, layer=args.layer)
    print("RankGauss transformation completed.")
    


print(adata.obs[args.obs_y].value_counts())
model=scKnockPath()
layer=args.layer if args.layer is not None else 'rankgauss' if args.rank_gauss else None
X,Xh,y,groups_list=model.prepare_data(adata=adata,obs_y=args.obs_y,
                                      layer=layer,genesets=args.geneset_path,
                                      gene_names=args.gene_names,gene_thresh=args.gene_thresh,class1=args.class1,class2=args.class2)

# generate knockoff data
file_name=f'{args.cell_type}_{args.class2}_vs_{args.class1}_knockoff_seed={args.knockoff_seed}.npz'

data_file = os.path.join(data_folder, file_name)
load=args.load_knockoff
np.random.seed(args.knockoff_seed)



if load and os.path.exists(data_file):
    print('loading data')
    data = np.load(data_file,allow_pickle=True)
    Xc = data['Xc']
    y = pd.Series(data['y'])
else:
    print('generating data')
    Xc,X,Xk = model.knockoff_sampler(X,Xh,method='sdp')
    if args.save_knockoff:
        np.savez(data_file, Xc=Xc, y=y)
        
sys.stdout.flush()
# alphas=np.logspace(-1,-3, num=20)
# y 有两类，都是str，将 y 转化为 0，1 数值类型
# y = y.astype(float)  # 确保 y 是数值类型以进行计算
print(np.unique(y, return_counts=True))  # 再次检查 y 的类别分布
y[y==args.class1] = 0
y[y==args.class2] = 1
y = y.astype(float)  # 转换为数值类型以进行计算
# 每次 cv 中模型分数如果下降超过 0.05 就停止训练,防止过拟合，可以设置上涨很少或者下降很多时停止
# model.fit_W(Xc, y, groups_list)  # 直接计算 W 统计量，
model.fit(Xc, y, groups_list,n_alphas=30,cv=3,l1_ratio=float(args.l1_ratio))

final_pathways=model.select_pathway(fdr=args.fdr,offset=0)

print(f'scKnockPath selects {len(final_pathways)} pathways')
print(final_pathways)

if args.save_model:
    if args.cell_type is not None:
        model.save_model(f'{data_folder}/scKnockPath_l1_ratio_{args.l1_ratio}_{args.class1}_{args.class2}_{args.obs_y}_{args.cell_type}_knockoff_seed={args.knockoff_seed}.pkl')
    else:
        model.save_model(f'{data_folder}/scKnockPath_l1_ratio_{args.l1_ratio}_{args.class1}_{args.class2}_{args.obs_y}_knockoff_seed={args.knockoff_seed}.pkl')