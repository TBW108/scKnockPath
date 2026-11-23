"""Small CUDA integration check; generates no distributable data."""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import json,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.special import expit
from sklearn.model_selection import KFold
root=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='sckp-smoke-') as tmp:
    tmp=Path(tmp);rng=np.random.default_rng(42);x=rng.normal(size=(90,24));y=np.tile([0,1],45).astype(str)
    a=ad.AnnData(x,obs=pd.DataFrame({'label':y},index=[f'c{i}' for i in range(90)]),var=pd.DataFrame(index=[f'g{i}' for i in range(24)]))
    a.write_h5ad(tmp/'input.h5ad');gmt=tmp/'pathways.gmt'
    gmt.write_text('\n'.join('\t'.join([f'P{i}','na']+[f'g{j}' for j in range(i*4,(i+1)*4)]) for i in range(6)))
    subprocess.run([sys.executable,str(root/'scripts/run_gpu.py'),'--mode','real','--adata',str(tmp/'input.h5ad'),'--gmt',str(gmt),'--label','label','--control','0','--case','1','--gene-thresh','1','--max-iter','60','--out',str(tmp/'out')],check=True,stdout=subprocess.DEVNULL)
    r=json.loads((tmp/'out/result.json').read_text());z=np.load(tmp/'out/coefficients.npz');b=z['coef'];groups=np.split(np.arange(24),[4,8,12,16,20]);w=np.array([abs(b[g]).mean()-abs(b[g+24]).mean() for g in groups]);table=pd.read_csv(tmp/'out/pathways.csv')
    assert np.allclose(w,table.W,atol=1e-14)
    flag=np.zeros(6,bool)
    for t in np.sort(np.abs(w[w!=0])):
        if (1+(w<=-t).sum())/max(1,(w>=t).sum())<=.2:flag=w>=t;break
    assert list(table.pathway[flag])==r['selected']
    # Reconstruct deterministic knockoffs and all held-out errors independently.
    sys.path.insert(0,str(root/'model'))
    import torch
    from model_gpu import scKnockPathGPU
    model=scKnockPathGPU(exact_prox=True);xx,xh,_,_=model.prepare_data(a,'label',str(gmt),gene_thresh=1,class1='0',class2='1')
    np.random.seed(42);design,_,_=model.knockoff_sampler(xx,xh)
    design=design.cpu().numpy();cv=np.load(tmp/'out/cv.npz');n_grid=len(cv['grid']);errors=np.zeros((n_grid,3))
    for f,(_,te) in enumerate(KFold(3,shuffle=True,random_state=42).split(y)):
        assert np.all(cv['fold'][te]==f)
        cols=np.arange(f,3*n_grid,3);prob=expit(design[te]@cv['coef'][:,cols]+cv['intercept'][cols]);errors[:,f]=np.mean((prob>=.5)!=y[te,None].astype(float),axis=0)
    assert np.array_equal(errors,cv['errors']);means=errors.mean(1)
    assert r['lambda_value']==max(cv['grid'][means==means.min()])
    assert r['cv_at_cap']==int((cv['n_iter']>=60).sum())
    counts=rng.poisson(3,size=(90,24)).astype(float)
    sim=ad.AnnData(counts,obs=pd.DataFrame({'cell_type':y},index=a.obs_names),var=a.var.copy())
    sim.uns['all_pathways']={f'P{i}':np.array([f'g{j}' for j in range(i*4,(i+1)*4)]) for i in range(6)}
    sim.uns['effect_pathways']=np.array(['P0']);sim.write_h5ad(tmp/'sim.h5ad')
    subprocess.run([sys.executable,str(root/'scripts/run_gpu.py'),'--mode','simulation','--adata',str(tmp/'sim.h5ad'),'--max-iter','60','--out',str(tmp/'simout')],check=True,stdout=subprocess.DEVNULL)
    sr=json.loads((tmp/'simout/result.json').read_text());scv=np.load(tmp/'simout/cv.npz')
    assert len(scv['grid'])==20 and sr['knockoff_seed']==88
    chosen=set(sr['selected']);assert sr['FDP']==len(chosen-{'P0'})/max(1,len(chosen)) and sr['Power']==len(chosen&{'P0'})

print('PASS: CUDA pipeline, deterministic CV errors, tie rule, coefficient/W/selection reconstruction, cap flags.')
