"""Summarize generated scKnockPath runs as mean +/- standard error."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--results',default='results/fig2');a=p.parse_args();rows=[]
for f in Path(a.results).glob('*/result.json'):
    r=json.loads(f.read_text());setting=f.parent.name.rsplit('_seed=',1)[0]
    rows.append(dict(setting=setting,seed=r['configuration']['seed'],FDP=r['FDP'],Power=r['Power'],n_selected=r['n_selected'],cv_at_cap=r['cv_at_cap'],final_at_cap=r['final_at_cap']))
if not rows:raise SystemExit('No completed results found.')
df=pd.DataFrame(rows);out=[]
for setting,g in df.groupby('setting'):
    row=dict(setting=setting,n=len(g))
    for metric,label in [('FDP','FDR'),('Power','Power')]:row[label]=g[metric].mean();row[label+'_SE']=g[metric].std(ddof=1)/np.sqrt(len(g))
    out.append(row)
df.to_csv(Path(a.results)/'per_seed.csv',index=False);pd.DataFrame(out).to_csv(Path(a.results)/'summary.csv',index=False);print(pd.DataFrame(out).to_string(index=False))
