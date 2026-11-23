# Run in an R >= 4.3 environment; scDesign3 installation may require system libraries.
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c("SingleCellExperiment", "DuoClustering2018", "zellkonverter", "scuttle"), ask=FALSE, update=FALSE)
install.packages(c("remotes", "tidyverse", "reticulate", "circlize", "GSA", "Matrix"), repos="https://cloud.r-project.org")
remotes::install_github("SONGDONGYUAN1994/scDesign3", upgrade="never")
