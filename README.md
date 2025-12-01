# scKnckPath
A statistical model to find differential pathways that are conditionally correlated with target phenotypes by deconfounding the effect of overlapping genes with FDR control. This framework requires the scRNA-seq data preprocessed by scanpy, target FDR and .gmt pathway database file as input and outputs a ranked list of pathways.

## 1. Create the Enviroment
`conda env create -f sckp_env.yml`

## 2. Preprocess the Dataset
For the quality-controlled gene count matrix, we sequentially applied 
`scanpy.pp.log1p()` and `scanpy.pp.normalize\_total()` to obtain 
normalized expression values. We then identified the top $G$ (2000 in default) highly variable genes 
using `scanpy.pp.highly\_variable\_genes()` and extracted their expression 
profiles.

## 2. Run the scKnockPath to Identify Pathways
Prepare the `.h5ad` scRNA-seq dataset and `.gmt` pathway database in the folder `./data`.
In the commandline, run the following example code to obtain the results:
```
python ./test/real_exp/scKnockPath_run.py \
            --adata ./data/T1D/processed_hpap.h5ad \
            --obs_y diabetes_status \
            --geneset_path ./data/pathway_db/reactome_kegg.gmt \
            --gene_thresh 15 \
            --class1 ND 
            --class2 AAB+
            --fdr 0.2 \
            --cell_type Beta
```

## 3. Output the Results
The model is saved as a `.pkl` file and you can upload it using 'pickle' package.
```
cell_type="Beta"
with open(f'results/real_exp/processed_hpap/scKnockPath_results/scKnockPath_ND_AAB+_{cell_type}.pkl', 'rb') as f:
        model = pickle.load(f)
selected_pathways=pd.DataFrame(model.rank_pathways()[::-1][:len(model.select_pathway(0.2))]
print(selected_pathways)
```

  
