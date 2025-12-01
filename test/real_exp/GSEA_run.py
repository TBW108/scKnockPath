import numpy as np
import pandas as pd
import gseapy as gp
import scanpy as sc
import os
import argparse
import sys


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
parser.add_argument('--cell_type', type=str,default=None, help='Cell type for analysis')


args = parser.parse_args()



def GSEA_test(adata,layer,gs_filepath,gene_names,obs_y,fdr,gene_thresh,save_path):

    adata.layers['lognorm'] = adata.X.layers[layer] if layer is not None else adata.X

    
    if hasattr(adata.layers['lognorm'], "toarray"):
        X = adata.layers['lognorm'].toarray()
    else:
        X = adata.layers['lognorm']
    # X = StandardScaler().fit_transform(X)
    
    if gene_names is not None: 
        adata.var_names = adata.var[gene_names].astype(str)
    
    df = pd.DataFrame(X, columns=adata.var_names)
    geneset_dict = gp.read_gmt(gs_filepath)
    
    res = gp.gsea(
        data=df.T,
        gene_sets=geneset_dict,
        cls=adata.obs[obs_y],
        permutation_num=1000,
        method='signal_to_noise',
        min_size=gene_thresh,
        seed=42
    )
    # print(res.res2d)
    result = res.res2d[(res.res2d['FDR q-val'] <= fdr)]['Term']
    
    if save_path is not None:
        pd.DataFrame(result, columns=['Term']).to_csv(save_path, index=False)
        print(f'GSEA results are saved to {save_path}')
    print(f'There are {len(result)} pathways selected with FDR<={fdr}')
    
    print('-------------------------------------------------------------------')
        
    
np.random.seed(42)

# Extract filename prefix from adata path
adata_filename = os.path.basename(args.adata)
filename_prefix = os.path.splitext(adata_filename)[0]

# Create corresponding folder in ./results
data_folder = f'./results/real_exp/{filename_prefix}'
os.makedirs(data_folder, exist_ok=True)
method_folder = os.path.join(data_folder, 'GSEA_results')
os.makedirs(method_folder, exist_ok=True)

adata=sc.read_h5ad(args.adata)
if args.cell_type is not None:
    if 'cell_type' not in adata.obs.columns:
        print("Error: 'cell_type' column not found in adata.obs. Please ensure your data contains cell type annotations.")
        sys.exit(1)
    adata = adata[adata.obs['cell_type'] == args.cell_type, :]
    print(f"Filtered data to {args.cell_type} cells: {adata.n_obs} cells remaining")

adata=adata[(adata.obs[args.obs_y]==args.class2)|(adata.obs[args.obs_y]==args.class1),:]
    
print(adata.obs[args.obs_y].value_counts())

save_path=f'./results/real_exp/{filename_prefix}/GSEA_results/GSEA_{args.class2}_vs_{args.class1}_{args.cell_type}.csv'

GSEA_test(adata=adata,layer=args.layer,gs_filepath=args.geneset_path,gene_names=args.gene_names, obs_y=args.obs_y,fdr=args.fdr, gene_thresh=args.gene_thresh,save_path=save_path)

