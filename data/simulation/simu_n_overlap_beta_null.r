basilisk_lib <- "/disk/tanbowen/.cache/R/basilisk/1.20.0/zellkonverter/1.18.0/zellkonverterAnnDataEnv-0.11.4/lib"

current_ld <- Sys.getenv("LD_LIBRARY_PATH")
ld_components <- c(basilisk_lib)
if (nzchar(current_ld)) {
  ld_components <- c(ld_components, strsplit(current_ld, ":")[[1]])
}
ld_components <- unique(ld_components)
Sys.setenv(LD_LIBRARY_PATH = paste(ld_components, collapse = ":"))

if (!any(grepl("libnghttp2.so", names(getLoadedDLLs())))) {
  dyn.load(file.path(basilisk_lib, "libnghttp2.so.14"), local = FALSE)
}

library(SingleCellExperiment)
library(tidyverse)
library(zellkonverter)
library(reticulate)
library(scuttle)
library(circlize)
set.seed(123)


# Example usage:
# seeds <- 46:55
seeds <- 46:75
n_overlaps <- seq(0, 10, 2)
for (n_overlap in n_overlaps) {
  dir.create("./null_overlap_data/", showWarnings = FALSE)
  path <- paste0("./null_overlap_data/n_overlap_beta", n_overlap)
  dir.create(path, showWarnings = FALSE)

  for (seed in seeds) {
    example_sce <- readH5AD(paste0("./overlap_data/n_overlap_beta", n_overlap, "/simu_scRNAseq_100pathways_", n_overlap, "overlaplog_seed=", seed, ".h5ad"))
    print(paste0("n_overlap=", n_overlap, ", seed=", seed))

    # extract the log-normalized assay
    log_example_newcount <- logcounts(logNormCounts(example_sce))

    # extract the DE pathways
    groups <- metadata(example_sce)$all_pathways
    DE_pathways <- metadata(example_sce)$effect_pathways
    print(paste("DE pathways:", paste(DE_pathways, collapse = ", ")))

    # generate a new y response with beta=0
    betas <- list()
    y <- 0
    for (pathway in DE_pathways) {
      effect_genes <- unlist(groups[pathway])
      zero_indices <- sample(seq_along(effect_genes), size = floor(length(effect_genes) * 0.2))

      betas[[pathway]] <- 0 * runif(length(effect_genes), -3, 3)
      betas[[pathway]][zero_indices] <- 0

      # x <- example_newcount[effect_genes, ]
      x <- log_example_newcount[effect_genes, ]
      y <- y + betas[[pathway]] %*% x
    }

    # Apply sigmoid function to y
    sigmoid <- function(x) {
      return(1 / (1 + exp(-x)))
    }
    y <- sigmoid(y - median(y) + rnorm(length(y), 0, 1))

    # save the new sce with null effect
    colData(example_sce)$cell_type <- as.vector(ifelse(y > 0.5, 1, 0))
    file_name <- paste0("./null_overlap_data/n_overlap_beta", n_overlap, "/simu_scRNAseq_100pathways_", n_overlap, "overlap_null_seed=", seed, ".h5ad")
    writeH5AD(example_sce, file_name)
    print(paste(file_name, "create success."))
  }
}
