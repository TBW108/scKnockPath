# Data sources and preparation

No expression matrices, generated datasets, pathway databases, or fitted models are distributed. Download data from their original providers and follow the providers' access terms. Links below were checked on 2026-11-21; portals may change export formats.

## Real expression data

- **HPAP:** [Gaulton Lab HPAP resource](https://www.gaultonlab.org/pages/Islet_expression_HPAP/) links expression downloads, donor annotations, and [analysis pipelines](https://github.com/Gaulton-Lab/HPAP-scRNA-seq). The current page offers an R object and annotations. The study preprocessing notebook expects the corresponding annotated `hpap_islet_scRNAseq.h5ad` export, with log-normalized `X` and `Cell Type`, `Diabetes Status`, `Sex` metadata. If downloading an R object, export it to AnnData while retaining that scale and metadata; an arbitrary raw-count export is not interchangeable. The release does not silently treat `expm1(log-normalized X)` as original counts.
- **HLCA:** [HLCA v1.0 download portal](https://data.humancellatlas.org/hca-bio-networks/lung/atlases/lung-v1-0), [source studies](https://data.humancellatlas.org/hca-bio-networks/lung/atlases/lung-v1-0/datasets), and [metadata documentation](https://github.com/LungCellAtlas/HLCA). Choose the **full/disease-inclusive** atlas, not the healthy core. The provider documents raw counts in `adata.raw.X`; the original `lung_disease.h5ad` workflow uses the annotated normalized expression in `X`. Required metadata: `assay`, `disease`, `cell_type`, `donor_id`; gene symbols in `var.feature_name`.
- **T2D:** [GEO GSE221156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156) supplies study data and sample metadata (17 ND, 17 T2D, 14 prediabetic donors). Use its corresponding annotated beta-cell count export, retaining `cell_type`, `disease`, `donor_id`, and `feature_name`. `real/prepare.py` expects this export as input; GEO raw matrices alone require metadata/annotation assembly first. The analysis uses all normal/T2D donors, not the historical female/age-restricted subset. The exact historical CELLxGENE export URL is not recorded in the source preprocessing code, so the study accession is supplied instead of an unverified asset URL.

## Preprocess

```bash
python real/prepare.py --dataset hpap --input data/hpap_islet_scRNAseq.h5ad --out data/hpap_processed.h5ad
python real/prepare.py --dataset hlca --input data/lung_disease.h5ad --out data/hlca_processed.h5ad
python real/prepare.py --dataset t2d --input data/beta_T2D.h5ad --out data/t2d_processed.h5ad
```

HPAP preserves the original normalized matrix and selects male Alpha/Beta cells from ND/AAB+ donors. HLCA selects 10x 3′ v3 normal/COVID-19 cells and 2,000 HVGs from the normalized matrix. T2D performs the original QC (≥200 genes/cell, genes in ≥3 cells, mitochondrial percentage <5), normalization to 10,000, log1p, 2,000 HVGs, and a seed-0 cap of 400 cells/donor after restricting to normal/T2D. Its mitochondrial rule uses `var_names` exactly as the original code; verify the export's identifier convention before interpreting that QC.

The output-free notebooks retain the original exploratory preprocessing history. For the revised T2D all-donor analysis use `prepare.py`, not the old age/sex restriction shown in `process.ipynb`. These commands require the stated annotated source exports; exact reproduction depends on matching their versions, gene identifiers, expression scale, and metadata.

The formal filtered inputs used for the main analyses had HPAP 15,235 cells/193 pathways, HLCA 13,602 cells/134 pathways, and T2D 13,600 cells/73 pathways (2,000 genes each). Inspect any discrepancy rather than silently changing filters to force these counts.

## Pathways and simulation reference

Download gene-symbol GMT files from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/). The study uses a combined human KEGG/Reactome library, filtered to at least 15 genes present in each real dataset. Save the same chosen release and record its checksum: changing database release changes candidate pathways. The original bundled GMT is intentionally not redistributed here.

Simulation scripts retrieve the reference through [DuoClustering2018](https://bioconductor.org/packages/DuoClustering2018/), then use [scDesign3](https://github.com/SONGDONGYUAN1994/scDesign3). See [SIMULATIONS.md](SIMULATIONS.md).
