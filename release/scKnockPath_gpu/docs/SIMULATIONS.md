# Simulation reproduction

Run commands from the release root. Install the Python/CUDA environment and R dependencies first. The reference `sce_filteredExpr10_Zhengmix4eq` is fetched by DuoClustering2018; its initial download needs network access. Record `sessionInfo()` and the scDesign3 commit for reproducibility: the installer follows upstream rather than claiming a historical frozen R environment.

The generators preserve the study's negative-binomial/Gaussian-copula construction: 1,500 genes partitioned into 80 pathways (40×15, 20×20, 20×25), 10 differential pathways, and DE genes drawn from each pathway's original members. The `own` scripts retain the original overlap-dependent copula adjustment and in-place neighbor-sharing loop. The overlap parameter is the number borrowed in each direction, **not necessarily the final pairwise intersection size**. These are the actual study generators, not a redesigned simulation.

## Generate the Fig. 2 inputs

```bash
Rscript --vanilla simulation/gen_overlap.r 0 46 95 4
Rscript --vanilla simulation/gen_overlap_own.r 5,10 46 95 4
Rscript --vanilla simulation/gen_signal_own.r 0.4,0.6,0.8 46 95 4
Rscript --vanilla simulation/gen_sample_own.r 600,1200,1800 46 95 4
Rscript --vanilla simulation/gen_sample_own.r 3000,4000,5000 46 65 4
```

Arguments are `values seed_from seed_to cores`. Set `SCKP_SIM_DIR=/path/to/simulation` to relocate outputs. Existing generated files are skipped. The historical filename contains `100pathways`, but the actual library contains **80** pathways; the name is preserved for compatibility. Different scenarios have independently generated inputs, including their default settings.

## Fit and summarize

```bash
python simulation/run_fig2.py --dry-run
CUDA_VISIBLE_DEVICES=0 python simulation/run_fig2.py --data-root data/simulation --out results/fig2
python simulation/summarize.py --results results/fig2
```

This schedules 510 scKnockPath fits: 50 seeds per standard setting, and 20 seeds for each of N=3,000/4,000/5,000. `--seed-from 46 --seed-to 46` is a small initial check. The batch runner skips completed outputs and stops on errors; a partially written run directory must be inspected before removing/retrying it.

The simulation entry point performs the same hk20 normalization and log1p transformation as the current GPU experiments. It uses the supplied `all_pathways` and `effect_pathways` annotations, knockoff seed `42 + data seed`, and KFold seed 42. It reports per-dataset FDP and power; the summary reports their means and SE (sample SD divided by √n). No extra `cv_gain` filtering is used. All CV endpoints and final coefficients are saved, including capped fits.

The embedded Fig. 2 is a rendering of the current manuscript figure and includes historical GSEA/SCPA/PADOG/CAMERA/SOGL/OGL results. This release's batch runner does not regenerate those comparator arms or claim their historical results use the revised scKnockPath solver. The figure's error bars are bootstrap 95% CIs, whereas `summarize.py` outputs SE.
