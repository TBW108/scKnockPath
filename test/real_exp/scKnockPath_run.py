import sys
import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
sys.path.append('./model/') # 添加model文件夹到系统路径
from util import pancreas_preprocess
from model_lambda import scKnockPath
import argparse
np.random.seed(42)

parser = argparse.ArgumentParser(description='scKnockPath pathway analysis')
parser.add_argument('--adata', type=str, required=True, help='Path to adata h5ad file')
parser.add_argument('--layer', type=str, default=None, help='Layer to use from adata')
parser.add_argument('--gene_names', type=str, default=None, help='Gene names obs (default: adata.var_names)')
parser.add_argument('--geneset_path', type=str, required=True, help='Path to pathway database file')
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

args = parser.parse_args()



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


print(adata.obs[args.obs_y].value_counts())
model=scKnockPath()
X,Xh,y,groups_list=model.prepare_data(adata=adata,obs_y=args.obs_y,
                                      layer=args.layer,genesets=args.geneset_path,
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


alphas=np.logspace(-1,-3, num=15)
# 直到cv score第一次下降时停止
model.fit(Xc,y,groups_list,alphas=alphas,thresh=0.01)
final_pathways=model.select_pathway(fdr=args.fdr,offset=0)


print(f'scKnockPath selects {len(final_pathways)} pathways')
print(final_pathways)

if args.save_model:
    if args.cell_type is not None:
        model.save_model(f'{data_folder}/scKnockPath_{args.class1}_{args.class2}_{args.cell_type}_knockoff_seed={args.knockoff_seed}.pkl')
    else:
        model.save_model(f'{data_folder}/scKnockPath_{args.class1}_{args.class2}_knockoff_seed={args.knockoff_seed}.pkl')