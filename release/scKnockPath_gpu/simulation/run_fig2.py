"""Run only the GPU scKnockPath arm of Fig2 on locally generated datasets."""
import argparse,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',default='data/simulation');p.add_argument('--out',default='results/fig2');p.add_argument('--seed-from',type=int,default=46);p.add_argument('--seed-to',type=int,default=95);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
root=Path(__file__).resolve().parents[1]
for scenario,values in [('overlap',[0,5,10]),('signal',[.4,.6,.8]),('sample',[600,1200,1800,3000,4000,5000])]:
    for value in values:
        native='overlap' if scenario=='overlap' and value==0 else scenario+'own'
        for seed in range(a.seed_from,min(a.seed_to,65) + 1 if scenario=='sample' and value>=3000 else a.seed_to+1):
            stem=f'{native}={value:g}';src=Path(a.data_root)/f'{native}_data'/stem/f'simu_scRNAseq_100pathways_{value:g}{native}_seed={seed}.h5ad'
            out=Path(a.out)/f'{scenario}={value:g}_seed={seed}'
            cmd=[sys.executable,str(root/'scripts/run_gpu.py'),'--mode','simulation','--adata',str(src),'--seed',str(seed),'--out',str(out)]
            if a.dry_run:print(' '.join(cmd));continue
            if (out/'result.json').exists():
                r=json.loads((out/'result.json').read_text())
                if r['configuration']['seed']!=seed or Path(r['configuration']['adata'])!=src:raise RuntimeError(f'Conflicting result: {out}')
                print('Already complete:',out);continue
            if not src.exists():raise FileNotFoundError(f'Generate this dataset first: {src}')
            subprocess.run(cmd,check=True)
