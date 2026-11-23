# GPU execution

Run from the release root after installing the environment. `run_gpu.py --help` lists options. Select a physical GPU with `CUDA_VISIBLE_DEVICES`; the process then addresses it as `cuda:0`. The solver uses float64. Memory usage grows with the number of cells, unfolded pathway memberships, and 60 simulation or 90 real-data simultaneous CV problems; no universal minimum GPU memory is claimed.

## Real analyses

After the preparation in [DATA.md](DATA.md), place the selected gene-symbol library at `data/pathways.gmt`:

```bash
CUDA_VISIBLE_DEVICES=0 python scripts/run_gpu.py --mode real \
  --adata data/hpap_processed.h5ad --gmt data/pathways.gmt \
  --label diabetes_status --control ND --case AAB+ --cell-type Beta \
  --out results/hpap

CUDA_VISIBLE_DEVICES=0 python scripts/run_gpu.py --mode real \
  --adata data/hlca_processed.h5ad --gmt data/pathways.gmt \
  --label disease --control normal --case COVID-19 \
  --cell-type 'elicited macrophage' --gene-names feature_name --out results/hlca

CUDA_VISIBLE_DEVICES=0 python scripts/run_gpu.py --mode real \
  --adata data/t2d_processed.h5ad --gmt data/pathways.gmt \
  --label disease --control normal --case 'type 2 diabetes mellitus' \
  --cell-type 'type B pancreatic cell' --gene-names feature_name --out results/t2d
```

The real-data entry point assumes `X` is already normalized/log1p and does not normalize it again. `--layer NAME` selects a different preprocessed expression layer. Real knockoff seed defaults to 42. Main analyses use random-cell KFold CV; `--donor-key donor_id` explicitly requests the separate donor-grouped sensitivity analysis and is not interchangeable with the main setting.

## Fixed fitting rule

The public entry point calls the GPU batched solver directly: exact proximal operator, α=0.3, 20 λ values for simulations or 30 for real data from 0.1 to 0.001 rounded to two significant digits, three-fold CV seed 42, classification at p≥0.5, minimum mean fold error, and largest λ for exact ties. Each CV/final fit stops at `--max-iter` (default 3,000), tol=1e-7; no automatic expansion and no `cv_gain` rule. Do not call the inherited historical `fit()` method to reproduce this entry point's rule.

Covariance options are `--covariance ledoitwolf` (default), `graphicallasso` (regularization 0.1), or `ridge --ridge-eps VALUE`. Ridge means sample covariance + εI; the GPU implementation computes the equivalent centered matrix product. The R3-M1 ε=1 experiment is an explicitly selected sensitivity configuration, not a changed model default.

## Outputs and interpretation

- `pathways.csv`: pathway names, original/knockoff scores, W, selection flag.
- `result.json`: configuration, input checksum, λ, selected pathways, convergence flags, construction diagnostics; simulation runs also have FDP and power.
- `coefficients.npz`: final coefficients/intercept and group sizes, sufficient to reconstruct W and selection.
- `cv.npz`: all CV coefficients/intercepts, iteration counts, errors, fold assignments, and λ grid.

The threshold is knockoff+ (offset 1), q=0.2 by default. An empty selected set has FDP=0. Power is undefined for global-null inputs. Preserve capped fits and report their flags. Finite-sample FDR statements require the knockoff assumptions; overlapping pathway experiments provide empirical evidence, not a blanket guarantee.

## Development check

`python tests/smoke_gpu.py` runs a small synthetic CUDA example in a temporary directory and independently reconstructs W, selection, and the CV choice. It is a software check, not a simulation benchmark. The release model files are unchanged copies; their source hashes are recorded in `docs/source_provenance.json`.
