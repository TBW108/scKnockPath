"""Portable CUDA entry point using the revised classification-error selection rule."""
import argparse
import os
from pathlib import Path
import sys
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'model'))

def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata',required=True);p.add_argument('--out',required=True)
    p.add_argument('--mode',choices=['simulation','real'],required=True)
    p.add_argument('--gmt');p.add_argument('--label',default='cell_type')
    p.add_argument('--control',default='0');p.add_argument('--case',default='1')
    p.add_argument('--gene-names');p.add_argument('--cell-type');p.add_argument('--cell-type-key',default='cell_type')
    p.add_argument('--layer');p.add_argument('--gene-thresh',type=int)
    p.add_argument('--seed',type=int,default=46,help='Simulation data seed; knockoff seed defaults to seed+42')
    p.add_argument('--knockoff-seed',type=int);p.add_argument('--cv-seed',type=int,default=42)
    p.add_argument('--donor-key',help='Optional donor-grouped CV; omit for the main random-cell analysis')
    p.add_argument('--l1-ratio',type=float,default=.3);p.add_argument('--q',type=float,default=.2)
    p.add_argument('--max-iter',type=int,default=3000);p.add_argument('--device',default='cuda:0')
    p.add_argument('--covariance',choices=['ledoitwolf','ridge','graphicallasso'],default='ledoitwolf')
    p.add_argument('--ridge-eps',type=float,default=.01)
    return p

def main():
    args=parser().parse_args()
    import json,hashlib
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import torch
    from sklearn.model_selection import KFold,StratifiedGroupKFold
    from model_gpu import scKnockPathGPU,BatchedLogisticSGL,_FittedSGL
    out=Path(args.out)
    if out.exists():raise SystemExit(f'Refusing to overwrite {out}; choose a new output directory.')
    if not torch.cuda.is_available():raise SystemExit('CUDA is unavailable. Install a CUDA-enabled PyTorch build and check the NVIDIA driver.')
    if not args.device.startswith('cuda'):raise SystemExit('This entry point is for CUDA devices.')
    if args.case==args.control or args.max_iter<1:raise SystemExit('Distinct labels and positive max-iter are required.')
    torch.cuda.set_device(args.device);torch.use_deterministic_algorithms(True)
    a=ad.read_h5ad(args.adata)
    if args.cell_type:a=a[a.obs[args.cell_type_key].astype(str)==args.cell_type].copy()
    a=a[a.obs[args.label].astype(str).isin([args.control,args.case])].copy()
    a.obs[args.label]=a.obs[args.label].astype(str)
    if set(a.obs[args.label])!={args.control,args.case}:raise ValueError('Both requested classes must be present.')
    if args.layer:a.X=a.layers[args.layer].copy()
    truth=None
    if args.mode=='simulation':
        counts=np.asarray(a.X.toarray() if hasattr(a.X,'toarray') else a.X,dtype=float)
        b=ad.AnnData(counts.copy());sc.pp.normalize_total(b);sc.pp.log1p(b)
        sc.pp.highly_variable_genes(b,flavor='seurat',n_top_genes=counts.shape[1])
        disp=np.nan_to_num(b.var.dispersions_norm.values,nan=-np.inf)
        totals=counts[:,disp<=np.quantile(disp,.8)].sum(1)
        if np.any(totals<=0):raise ValueError('Zero hk20 library size; inspect generated counts.')
        a.X=np.log1p(counts/(totals/np.median(totals))[:,None])
        pathways={k:np.asarray(v) for k,v in a.uns['all_pathways'].items()}
        truth=set(np.asarray(a.uns['effect_pathways']).reshape(-1).astype(str))
    else:
        if not args.gmt:raise ValueError('--gmt is required for real data; input X must already be normalized/log1p.')
        pathways=args.gmt
    model=scKnockPathGPU(device=args.device,exact_prox=True)
    x,xh,labels,groups=model.prepare_data(a,obs_y=args.label,genesets=pathways,gene_names=args.gene_names,
        gene_thresh=args.gene_thresh or (1 if args.mode=='simulation' else 15),class1=args.control,class2=args.case)
    if len(groups)==0:raise ValueError('No pathway passes the gene threshold.')
    y=(labels==args.case).astype(float)
    kseed=args.knockoff_seed if args.knockoff_seed is not None else (args.seed+42 if args.mode=='simulation' else 42)
    np.random.seed(kseed)
    design,xs,xk=model.knockoff_sampler(x,xh,cov_estimator=args.covariance,ridge=args.ridge_eps)
    n_grid=20 if args.mode=='simulation' else 30
    grid=np.array([float(f'{v:.2g}') for v in np.logspace(-1,-3,n_grid)])
    if args.donor_key:
        splits=list(StratifiedGroupKFold(3,shuffle=True,random_state=args.cv_seed).split(y,y,a.obs[args.donor_key].astype(str)))
    else:splits=list(KFold(3,shuffle=True,random_state=args.cv_seed).split(y))
    d=sum(map(len,groups));solver=BatchedLogisticSGL(design,list(groups)+[g+d for g in groups],args.l1_ratio,exact_prox=True)
    target=torch.as_tensor(y,device=args.device,dtype=torch.float64)
    weights=torch.zeros((len(y),3*n_grid),device=args.device,dtype=torch.float64);fold=np.empty(len(y),int)
    for f,(tr,te) in enumerate(splits):
        if len(np.unique(y[tr]))!=2:raise ValueError('A training fold has only one class.')
        fold[te]=f
        for i in range(n_grid):weights[torch.as_tensor(tr,device=args.device),3*i+f]=1/len(tr)
    b,ic,nit=solver.solve(target,weights,np.repeat(grid,3),max_iter=args.max_iter,tol=1e-7)
    errors=np.empty((n_grid,3))
    for f,(_,te) in enumerate(splits):
        cols=np.arange(f,3*n_grid,3)
        probs=torch.sigmoid(design[torch.as_tensor(te,device=args.device)]@b[:,cols]+ic[cols]).cpu().numpy()
        errors[:,f]=np.mean((probs>=.5)!=y[te,None],axis=0)
    means=errors.mean(1);best=float(max(grid[means==means.min()]))
    out.mkdir(parents=True)
    np.savez_compressed(out/'cv.npz',coef=b.cpu().numpy(),intercept=ic.cpu().numpy(),n_iter=nit.cpu().numpy(),errors=errors,fold=fold,grid=grid)
    cvnit=nit.cpu().numpy();del b,ic,nit,weights
    b,ic,nit=solver.solve(target,torch.full((len(y),1),1/len(y),device=args.device,dtype=torch.float64),np.array([best]),max_iter=args.max_iter,tol=1e-7)
    coef=b[:,0].cpu().numpy();intercept=float(ic[0]);finalnit=int(nit[0])
    model.unfold_length=d
    model.update_Wg(_FittedSGL(coef,intercept,np.array([0.,1.])));chosen=model.select_pathway(args.q,offset=1,stat='W').tolist()
    pd.DataFrame(dict(pathway=model.pathway_names,Z=model.Z,Zk=model.Zk,W=model.Wg,selected=np.isin(model.pathway_names,chosen))).to_csv(out/'pathways.csv',index=False)
    np.savez_compressed(out/'coefficients.npz',coef=coef,intercept=intercept,sizes=np.array([len(g) for g in groups]))
    result=dict(configuration=vars(args),knockoff_seed=kseed,lambda_value=best,cv_mean_errors=means.tolist(),selected=chosen,n_selected=len(chosen),
        cv_at_cap=int((cvnit>=args.max_iter).sum()),final_at_cap=finalnit>=args.max_iter,final_n_iter=finalnit,
        q=args.q,offset=1,rule='minimum mean classification error, p>=0.5; exact ties choose largest lambda',knockoff_info=model.knockoff_info,
        input_sha256=hashlib.sha256(Path(args.adata).read_bytes()).hexdigest(),n_cells=len(y),n_genes=a.n_vars)
    if truth is not None:result.update(FDP=len(set(chosen)-truth)/max(1,len(chosen)),Power=len(set(chosen)&truth)/len(truth) if truth else None)
    (out/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
