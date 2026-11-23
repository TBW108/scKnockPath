"""Dataset-specific preprocessing; choose the correct upstream expression scale."""
import argparse
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset',choices=['t2d','hpap','hlca'],required=True);p.add_argument('--input',required=True);p.add_argument('--out',required=True);a=p.parse_args()
import scanpy as sc
import numpy as np
x=sc.read_h5ad(a.input)
if Path(a.out).exists():raise SystemExit('Output exists; refusing to overwrite.')
if a.dataset=='t2d':
    # Input must be the beta-cell raw-count export, as in make_alldonors.py.
    x=x[x.obs.cell_type=='type B pancreatic cell'].copy()
    x.var['mt']=x.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(x,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
    sc.pp.filter_cells(x,min_genes=200);sc.pp.filter_genes(x,min_cells=3)
    x=x[x.obs.pct_counts_mt<5].copy();sc.pp.normalize_total(x,target_sum=1e4);sc.pp.log1p(x)
    sc.pp.highly_variable_genes(x,n_top_genes=2000);x=x[:,x.var.highly_variable].copy()
    x=x[x.obs.disease.isin(['normal','type 2 diabetes mellitus'])].copy();rng=np.random.default_rng(0);keep=[]
    for _,idx in x.obs.groupby('donor_id',observed=True).indices.items():keep.extend(rng.choice(idx,size=min(len(idx),400),replace=False))
    x=x[np.sort(keep)].copy()
elif a.dataset=='hpap':
    # Upstream X is already log-normalized; never normalize it again.
    x=x[x.obs['Cell Type'].isin(['Alpha','Beta']) & x.obs['Diabetes Status'].isin(['ND','AAB+'])].copy()
    x.obs['cell_type']=x.obs.pop('Cell Type');x.obs['diabetes_status']=x.obs.pop('Diabetes Status')
    x=x[x.obs.Sex=='M'].copy()
else:
    # Matches the first saved HLCA branch; X is the upstream normalized expression.
    x=x[(x.obs.assay=="10x 3' v3") & x.obs.disease.isin(['normal','COVID-19'])].copy()
    sc.pp.highly_variable_genes(x,n_top_genes=2000,subset=True)
Path(a.out).parent.mkdir(parents=True,exist_ok=True);x.write_h5ad(a.out);print(x)
