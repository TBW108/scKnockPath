library("SeuratDisk")
# 读取Seurat文件
data_R <- readRDS("data/T1D/hpap_islet_scRNAseq.rds")

SaveH5Seurat(data_R, filename = "data/T1D/hpap_islet_scRNAseq.h5Seurat")

# 在目录获得h5ad文件
Convert("data/T1D/hpap_islet_scRNAseq.h5Seurat", dest = "h5ad")
