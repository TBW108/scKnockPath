library(scDesign3, verbose = FALSE)
library(SingleCellExperiment)
library(DuoClustering2018)
library(tidyverse)
library(zellkonverter)
library(reticulate)
library(scuttle)
library(circlize)
library(GSA)
set.seed(123)

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



# get data
Zhengmix4eq_sce <- get("sce_filteredExpr10_Zhengmix4eq")(metadata = FALSE)

# replace the row names with gene symbols
rownames(Zhengmix4eq_sce) <- rowData(Zhengmix4eq_sce)$symbol
example_sce <- Zhengmix4eq_sce

# Filter out lowly expressed genes and cells
example_sce <- example_sce[!is.na(rownames(example_sce)), ]
example_sce <- example_sce[rowSums(counts(example_sce) > 0) >= 5, ]
example_sce <- example_sce[, colSums(counts(example_sce)) > 100]

example_sce

# example_simu <- scdesign3(
#     sce = example_sce,
#     assay_use = "counts",
#     celltype = "phenoid",
#     pseudotime = NULL,
#     spatial = NULL,
#     other_covariates = NULL,
#     mu_formula = "phenoid",
#     sigma_formula = "1", # If you want your dispersion also varies along pseudotime, use "s(pseudotime, k = 5, bs = 'cr')"
#     family_use = "nb",
#     n_cores = 100,
#     usebam = FALSE,
#     corr_formula = "1",
#     copula = "gaussian",
#     DT = TRUE,
#     pseudo_obs = FALSE,
#     return_model = FALSE,
#     nonzerovar = FALSE,
#     ncell = 8000
#   )

# # extract T-cells
# selected_cells <- which(example_simu$new_covariate[,1] == c("regulatory.t"))
# example_newcount <- example_simu$new_count[,selected_cells]
# dim(example_newcount)

# names(example_simu)
# table(example_simu$new_covariate[,1])


# # for n_overlap in 0:20 step=2
# n_overlap=5
# # divide the genes into 100 pathways
# all_genes <- rowData(example_sce)$symbol
# num_pathways <- 100
# genes_per_pathway <- floor(length(all_genes) / num_pathways)
# groups <- split(all_genes, rep(1:num_pathways, each = genes_per_pathway, length.out = length(all_genes)))
# names(groups) <- paste0("pathway", seq_along(groups))
# groups

# # add n_overlap genes to each pathway
# if (n_overlap!=0)
#   {for (i in seq_along(groups)) {
#       if (i != length(groups)) {
#         adjacent_genes <- groups[[i + 1]][1:n_overlap]
#         groups[[i]] <- unique(c(groups[[i]], adjacent_genes))
#         }
#       }
#   }

# table(sapply(groups,length))
# seed=42
# # select 10 pathways as DE pathways
# set.seed(seed)
# DE_pathways <- sample(names(groups), 10)
# print(DE_pathways)

# # generate beta coefficients(with +- ) for 80% DE genes in DE pathways
# betas <- list()
# y=0
# for (pathway in DE_pathways) {
#     effect_genes = unlist(groups[pathway])
#     set.seed(seed)
#     zero_indices <- sample(seq_along(effect_genes), size = floor(length(effect_genes) * 0.2))
#     betas[[pathway]] <- runif(length(effect_genes), -3, 3)
#     betas[[pathway]][zero_indices] <- 0

#     # generate epsilon for each gene
#     # epsilons[[pathway]] <- rnorm(length(effect_genes), 0, 1)

#     x = example_newcount[effect_genes,]
#     y <- y + betas[[pathway]]%*% x
#     dim(y)
# }

# dim(y)
# # Calculate the mean and median of y
# y_mean <- mean(y)
# y_median <- median(y)

# print(paste("Mean of y:", y_mean))
# print(paste("Median of y:", y_median))

# # apply sigmoid function to y
# sigmoid <- function(x) {
#   return(1 / (1 + exp(-x)))
# }
# y=sigmoid(y-y_median+rnorm(length(y), 0, 1))
# # rnorm(length(y), 0, 1)

# # Convert y to binary data
# threshold <- 0.5
# y_binary <- ifelse(y > threshold, 1, 0)

# new_sce <- SingleCellExperiment(assays = list(counts = example_newcount),
#                                   colData = list(cell_type = as.vector(y_binary)),
#                                   rowData = list(gene_name = rownames(example_sce)),
#                                   metadata = list(effect_pathways = DE_pathways, all_pathways = groups,true_betas=betas))

# colData(new_sce)
# # table(y_binary)
# # colData(example_sce)$cell_type <- as.factor(y_binary)
# # metadata(example_sce)$all_pathways <- groups
# # metadata(example_sce)$DE_pathways <- DE_pathways
# # Save the example_sce object as an h5ad file
# path=paste0("/data/tanbowen/scKnockPath/simulation_exp/gen_data/n_overlap_beta",n_overlap)
# dir.create(path, showWarnings = FALSE)
# file_name <- paste0(path,'/simu_scRNAseq_100pathways_', n_overlap, 'overlap_seed=', seed, '.h5ad')
# writeH5AD(new_sce, file_name)


generate_sce_with_sparsity <- function(example_sce, n_sparsity, DE_pathways, groups, seed, n_overlap = 4, save = TRUE) {
  set.seed(seed)
  example_simu <- scdesign3(
    sce = example_sce,
    assay_use = "counts",
    celltype = "phenoid",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    mu_formula = "phenoid",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 100,
    usebam = FALSE,
    corr_formula = "1",
    copula = "gaussian",
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = FALSE,
    ncell = 8000
  )

  # extract T-cells
  selected_cells <- which(example_simu$new_covariate[, 1] == c("regulatory.t"))
  example_newcount <- example_simu$new_count[, selected_cells]
  print(dim(example_newcount))
  print("------------------------------------------------")


  # Add n_overlap genes to each pathway
  if (n_overlap != 0) {
    for (i in seq_along(groups)) {
      if (i != length(groups)) {
        adjacent_genes <- groups[[i + 1]][1:n_overlap]
        groups[[i]] <- unique(c(groups[[i]], adjacent_genes))
      }
    }
  }

  # Generate beta coefficients for 80% DE genes in DE pathways
  betas <- list()
  y <- 0
  temp_sce <- SingleCellExperiment(assays = list(counts = example_newcount))
  log_example_newcount <- logcounts(logNormCounts(temp_sce))

  print(log_example_newcount[1:10, 1:10])
  for (pathway in DE_pathways) {
    effect_genes <- unlist(groups[pathway])
    zero_indices <- sample(seq_along(effect_genes), size = floor(length(effect_genes) * 0.2))
    betas[[pathway]] <- runif(length(effect_genes), -3, 3)
    betas[[pathway]][zero_indices] <- 0

    x <- log_example_newcount[effect_genes, ]
    y <- y + betas[[pathway]] %*% x
  }

  # Apply sigmoid function to y
  sigmoid <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  y <- sigmoid(y - median(y) + rnorm(length(y), 0, 1))

  # Convert y to binary data
  threshold <- 0.5
  y_binary <- ifelse(y > threshold, 1, 0)

  if (save) {
    new_sce <- SingleCellExperiment(
      assays = list(counts = example_newcount),
      colData = list(cell_type = as.vector(y_binary)),
      rowData = list(gene_name = rownames(example_sce)),
      metadata = list(effect_pathways = DE_pathways, all_pathways = groups, true_betas = betas)
    )
    # path <- paste0("/data/tanbowen/scKnockPath/simulation_exp/gen_data/n_overlap_beta", n_overlap)
    path <- paste0("/disk/tanbowen/scKnockPath/data/simulation/sparsity_data/n_sparsity_beta", length(DE_pathways))
    dir.create(path, showWarnings = FALSE)
    file_name <- paste0(path, "/simu_scRNAseq_100pathways_", length(DE_pathways), "sparsity_seed=", seed, ".h5ad")
    writeH5AD(new_sce, file_name)
    print(paste(file_name, "create success."))
  }
}

set.seed(42)
# Divide the genes into 100 pathways
all_genes <- rowData(example_sce)$symbol
num_pathways <- 100
genes_per_pathway <- floor(length(all_genes) / num_pathways)
groups <- split(all_genes, rep(1:num_pathways, each = genes_per_pathway, length.out = length(all_genes)))
names(groups) <- paste0("pathway", seq_along(groups))

# Select 10 pathways as DE pathways
# DE_pathways <- sample(names(groups), 10)
# print(DE_pathways)

# Example usage:
# seeds <- 46:55
seeds <- 56:75
n_sparsitys <- seq(6, 16, 2)
for (n_sparsity in n_sparsitys) {
  set.seed(42)
  DE_pathways <- sample(names(groups), n_sparsity)
  print(DE_pathways)
  for (seed in seeds) {
    generate_sce_with_sparsity(example_sce, n_sparsity, DE_pathways, groups, seed, n_overlap = 4, save = TRUE)
  }
}







# generate cell_type by using the formula cell_type=sigmoid(beta*x+epsilon)

# generate new data by cell_type

# save the data


# # get parameters
# example_data <- construct_data(
#   sce = example_sce,
#   assay_use = "counts",
#   celltype = "cell_type",
#   pseudotime = NULL,
#   spatial = NULL,
#   other_covariates = NULL,
#   corr_by = "1",
#   ncell = 2000,
# )

# example_marginal <- fit_marginal(
#   data = example_data,
#   predictor = "gene",
#   mu_formula = "cell_type",
#   sigma_formula = "1",
#   family_use = "nb",
#   n_cores = 100,
#   usebam = FALSE
# )

# set.seed(123)
# example_copula <- fit_copula(
#   sce = example_sce,
#   assay_use = "counts",
#   marginal_list = example_marginal,
#   family_use = "nb",
#   copula = "gaussian",
#   n_cores = 100,
#   input_data = example_data$dat,
# )

# example_para <- extract_para(
#   sce = example_sce,
#   assay_use = "counts",
#   marginal_list = example_marginal,
#   n_cores = 100,
#   family_use = "nb",
#   new_covariate = example_data$newCovariate,
#   data = example_data$dat
# )


# example_newcount <- simu_new(
#     sce = example_sce,
#     mean_mat = example_para$mean_mat,
#     sigma_mat = example_para$sigma_mat,
#     zero_mat = example_para$zero_mat,
#     quantile_mat = NULL,
#     copula_list = example_copula$copula_list,
#     n_cores = 100,
#     family_use = "nb",
#     input_data = example_data$dat,
#     new_covariate = example_data$newCovariate,
#     important_feature = example_copula$important_feature,
#     filtered_gene = example_data$filtered_gene
#   )

# table(sapply(groups,length))

# new_sce <- SingleCellExperiment(assays = list(counts = example_newcount),
#                                   colData = list(cell_type = as.vector(example_data$newCovariate$cell_type)),
#                                   rowData = list(gene_name = rownames(example_sce)),
#                                   metadata = list(DE_pathways = DE_pathways, all_pathways = groups,betas=betas))

# path=paste0("/data/tanbowen/scKnockPath/simulation_exp/gen_data/n_overlap_beta",n_overlap)
#     dir.create(path, showWarnings = FALSE)
# file_name <- paste0(path,'/simu_scRNAseq_100pathways_', n_overlap, 'overlap_seed=', seed, '.h5ad')
#     # file_name <- paste0(path,'/simu_scRNAseq_100pathways_', n_overlap, 'random_overlap_seed=', seed, '.h5ad')

#     writeH5AD(new_sce, file_name)
#     print(paste(file_name, 'create success.'))






# # Perform t-test for one of the effect_genes
# effect_gene <- effect_genes[10]
# effect_gene<- rownames(example_sce)[1800]
# gene_expression <- example_newcount[effect_gene, ]

# # Perform t-test
# cell_type <- example_data$newCovariate$cell_type
# t_test_result <- t.test(gene_expression ~ cell_type)
# print(t_test_result)




# generate_sce_with_overlap <- function(n_overlap,groups,pathways,seed=42,save=TRUE) {

#   # add overlapping genes
#   if (n_overlap!=0)
#   {for (i in seq_along(groups)) {
#       if (i != length(groups)) {
#         adjacent_genes <- groups[[i + 1]][1:n_overlap]
#         groups[[i]] <- unique(c(groups[[i]], adjacent_genes))
#         }
#       }
#   }

#   # for (i in seq_along(groups)) {
#   #   other_groups <- setdiff(seq_along(groups), i)
#   #   overlap_genes <- unlist(lapply(other_groups, function(j) sample(groups[[j]], n_overlap, replace = FALSE)))
#   #   groups[[i]] <- unique(c(groups[[i]], overlap_genes))
#   # }

#   ground_truth <- rep(FALSE, length(groups))
#   gene_exp_df <- tibble(gene = all_genes, expressed = rep(FALSE, length(all_genes)),DE_pathways=rep('',length(all_genes)))
#   new_cell_type <- unlist(as.vector(example_data$newCovariate['cell_type']))

#   for (pathway in pathways) {
#     ground_truth[which(names(groups) == pathway)] = TRUE
#     effect_genes = groups[pathway]
#     n <- length(unlist(effect_genes))
#     sample_size <- ceiling(n * 0.8)
#     set.seed(seed)
#     effect_genes <- sample(unlist(effect_genes), sample_size)

#     for (gene in unlist(effect_genes)) {
#       state <- sample(c(-1, 1), 1)
#       base <- mean(example_para$mean_mat[(new_cell_type == 0), gene])

#       if (gene_exp_df[gene_exp_df$gene == gene, 2] == FALSE) {
#         gene_exp_df[gene_exp_df$gene == gene, 2] = TRUE
#         gene_exp_df[gene_exp_df$gene == gene, 3] = pathway
#         if (state == 1) {
#           example_para$mean_mat[(new_cell_type == 1), gene] <- base * runif(1, 1.5, 2)
#         } else {
#           example_para$mean_mat[(new_cell_type == 1), gene] <- base * runif(1, 1/2, 2/3)
#         }
#       }
#       else{

#       }
#     }
#   }

#   for (gene in gene_exp_df$gene[!gene_exp_df$expressed]) {
#     mean_expr <- mean(example_para$mean_mat[, gene])
#     example_para$mean_mat[, gene] <- mean_expr
#   }

#   set.seed(seed)
#   example_newcount <- simu_new(
#     sce = example_sce,
#     mean_mat = example_para$mean_mat,
#     sigma_mat = example_para$sigma_mat,
#     zero_mat = example_para$zero_mat,
#     quantile_mat = NULL,
#     copula_list = example_copula$copula_list,
#     n_cores = 100,
#     family_use = "nb",
#     input_data = example_data$dat,
#     new_covariate = example_data$newCovariate,
#     important_feature = example_copula$important_feature,
#     filtered_gene = example_data$filtered_gene
#   )
#   print(example_newcount[1:10, 1:10])
#   print(table(gene_exp_df$expressed))
#   print(table(gene_exp_df$DE_pathways))
#   # Find genes with DE_pathways equal to "pathway100"
#   pathway49_genes <- gene_exp_df %>% filter(DE_pathways == "pathway49")
#   print(pathway49_genes)
#   if (save) {
#     new_sce <- SingleCellExperiment(assays = list(counts = example_newcount),
#                                   colData = list(cell_type = as.vector(new_cell_type)),
#                                   rowData = list(gene_name = rownames(example_sce)),
#                                   metadata = list(effect_pathways = pathways, all_pathways = groups,DE_genes=gene_exp_df))
#     path=paste0("/data/tanbowen/scKnockPath/simulation_exp/gen_data/n_overlap",n_overlap)
#     dir.create(path, showWarnings = FALSE)
#     file_name <- paste0(path,'/simu_scRNAseq_100pathways_', n_overlap, 'overlap_seed=', seed, '.h5ad')
#     # file_name <- paste0(path,'/simu_scRNAseq_100pathways_', n_overlap, 'random_overlap_seed=', seed, '.h5ad')

#     writeH5AD(new_sce, file_name)
#     print(paste(file_name, 'create success.'))
#   }

# }


# all_genes <- rowData(example_sce)$symbol
# num_pathways <- 100
# genes_per_pathway <- floor(length(all_genes) / num_pathways)
# groups <- split(all_genes, rep(1:num_pathways, each = genes_per_pathway, length.out = length(all_genes)))
# names(groups) <- paste0("pathway", seq_along(groups))
# set.seed(42)
# DE_pathways <- sample(names(groups), 10)
# print(DE_pathways)


# orig_mean_matrix<-example_para$mean_mat
seeds <- 46:55
# seeds=56:65
# Example usage:
for (n_overlap in 0:10) {
  for (seed in seeds) {
    example_para$mean_mat <- orig_mean_matrix
    generate_sce_with_overlap(n_overlap, groups, DE_pathways, seed = seed, save = TRUE)
  }
}


# # Example usage:
# orig_mean_matrix<-example_para$mean_mat
# for (n_overlap in 6:10) {
#   example_para$mean_mat=orig_mean_matrix
#   generate_sce_with_overlap(n_overlap)
# }

# n_overlap=5
# names(groups) <- paste0("pathway", seq_along(groups))
# for (i in seq_along(groups)) {
#   if (i !=length(groups)) {
#     adjacent_genes <- groups[[i + 1]][1:n_overlap]
#   } else{
#     adjacent_genes <- groups[[i - 1]][(length(groups[[i - 1]]) - n_overlap + 1):length(groups[[i - 1]])]
#   }
#   groups[[i]] <- unique(c(groups[[i]], adjacent_genes))
# }
# # groups
# pathway_num=sapply(groups,length)


# set.seed(123)
# pathways <- sample(names(groups), 10)

# ground_truth=rep(FALSE,length(groups))# labels of pathways

# # Create a data frame with gene names and their expression status
# gene_exp_df <- tibble(
#   gene = all_genes,
#   expressed = rep(FALSE, length(all_genes))
# )

# gene_exp_df
# table(gene_exp_df$expressed)
# new_cell_type=unlist(as.vector(example_data$newCovariate['cell_type']))

# for(pathway in pathways){
#   ground_truth[which(names(groups)==pathway)]=TRUE
#   effect_genes=groups[pathway]

#   # get 90% modified genes
#   n <- length(unlist(effect_genes))
#   sample_size <- ceiling(n * 0.9)
#   effect_genes <- sample(unlist(effect_genes), sample_size)

#   for(gene in unlist(effect_genes)){

#     state=sample(c(-1,1),1) # up or down regulation
#     base=mean(example_para$mean_mat[(new_cell_type==0), gene])

#     if(gene_exp_df[gene_exp_df$gene == gene,2]==FALSE){
#       gene_exp_df[gene_exp_df$gene == gene,2]=TRUE
#       if(state==1){
#         example_para$mean_mat[(new_cell_type==1), gene] <- base*runif(1,1.5,2)
#       }
#       else{
#         example_para$mean_mat[(new_cell_type==1), gene] <- base*runif(1,1/2,2/3)
#       }

#     }
#   }
# }

# # average the expression of non-DE genes
# for (gene in gene_exp_df$gene[!gene_exp_df$expressed]) {
#   mean_expr <- mean(example_para$mean_mat[, gene])
#   example_para$mean_mat[, gene] <- mean_expr
# }

# # table(gene_exp_df$expressed)
# # example_para$mean_mat[990:1015,effect_genes]
# # table(new_cell_type)

# # Check if example_para$mean_mat has any zero values
# # any(example_para$mean_mat <0)
# example_newcount <- simu_new(
#   sce = example_sce,
#   mean_mat = example_para$mean_mat,
#   sigma_mat = example_para$sigma_mat,
#   zero_mat = example_para$zero_mat,
#   quantile_mat = NULL,
#   copula_list = example_copula$copula_list,
#   n_cores = 50,
#   family_use = "nb",
#   input_data = example_data$dat,
#   new_covariate = example_data$newCovariate,
#   important_feature = example_copula$important_feature,
#   filtered_gene = example_data$filtered_gene
# )

# new_sce=SingleCellExperiment(assays=list(counts=example_newcount),
#                       colData=list(cell_type=as.vector(new_cell_type)),
#                       rowData=list(gene_name=rownames(example_sce)),
#                       metadata=list(effect_pathways=pathways,
#                                     all_pathways=groups)
#           )

# metadata(new_sce)$effect_pathways
# # metadata(new_sce)$effect_pathways
# # save results
# # file_name=paste0('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_10pathways_geneth=15.h5ad')

# # file_name=paste0('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_10pathways_geneth=20.h5ad') # nolint

# # file_name=paste0('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_10pathways_overlapadd.h5ad')
# # file_name=paste0('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_5pathways_overlapadd.h5ad')

# # file_name=paste0('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_10pathways_overlapadd100.h5ad')

# file_name=paste0('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_100pathways_3overlap.h5ad')

# # file_name=paste0('/data/tanbowen/scKnockPath/simulation_exp/gen_data/simu_scRNAseq_10pathways.h5ad')
# writeH5AD(new_sce, file_name)
# print(paste(file_name,'create success.'))

# test_sce=readH5AD(file_name)

# test_sce


# The first 5 groups are effective
# set.seed(123)
# for(n_effect_pathway in 1:5){
#     # select pathways
#     ground_truth=rep(FALSE,length(groups))
#     ground_truth[sample(1:length(groups), n_effect_pathway, replace = FALSE)]=TRUE
#     effect_pathways=groups[ground_truth]
#     # effect_pathways

#     # change in pathway level in mean matrix
#     for(pathway_genes in effect_pathways){
#       example_para$mean_mat[(new_cell_type==1), pathway_genes]<-mean(example_para$mean_mat[(new_cell_type==0),pathway_genes])*runif(1,4,5)
#     }

#     # generate new data # nolint: indentation_linter.
#     set.seed(123)
#     example_newcount <- simu_new(
#       sce = example_sce,
#       mean_mat = example_para$mean_mat,
#       sigma_mat = example_para$sigma_mat,
#       zero_mat = example_para$zero_mat,
#       quantile_mat = NULL,
#       copula_list = example_copula$copula_list,
#       n_cores = 50,
#       family_use = "nb",
#       input_data = example_data$dat,
#       new_covariate = example_data$newCovariate,
#       important_feature = example_copula$important_feature,
#       filtered_gene = example_data$filtered_gene
#     )
#     new_sce=SingleCellExperiment(assays=list(counts=example_newcount),
#                             colData=list(cell_type=as.vector(new_cell_type)),
#                             rowData=list(gene_name=rownames(example_sce)),
#                             metadata=list(pathways=groups,
#                                           ground_truth=ground_truth,
#                                           effect_pathways=names(effect_pathways)))

#     # save results
#     file_name=paste0('scKnockPath/compare_results/simulation/simu_scRNAseq',n_effect_pathway,'.h5ad')
#     for(t in 1:2){
#         writeH5AD(new_sce, file_name)
#     }
#     print(paste0(file_name,'create success.'))
# }

# "pathway31" "pathway15" "pathway14" "pathway3"  "pathway42"

# effect_genes=effect_pathways%>%
#     unlist()%>%
#     unique()



# change in single gene level
# for(gene in colnames(example_para$mean_mat)){
#   if(gene %in% effect_genes){
#     example_para$mean_mat[(new_cell_type==1), gene]<-mean(example_para$mean_mat[(new_cell_type==0), gene])*runif(1,5,10)
#   }
#   else{
#     example_para$mean_mat[(new_cell_type==1), gene]<-mean(example_para$mean_mat[(new_cell_type==0), gene])*runif(1,0.95,1.05)
#   }
# }

# example_para$mean_mat[1520:1530,effect_genes]

# effect_gene_id=rowData(example_sce)$symbol %in% effect_genes
# noneffect_gene_id=!effect_gene_id

# fold_change=2
# example_para$mean_mat[,effect_gene_id] <- apply(example_para$mean_mat[,effect_gene_id], 2, function(x){
#      max_mean=max(x,na.rm = TRUE)
#      x[x==max_mean]=x[x==max_mean]+effect
#      return(x)})

# example_para$mean_mat[,noneffect_gene_id] <- apply(example_para$mean_mat[,noneffect_gene_id], 2, function(x){
#     avg <- (max(x,na.rm = TRUE)+min(x,na.rm = TRUE))/2
#     new_mean <- rep(avg, length(x))
#      return(new_mean)})

# colMeans(example_para$mean_mat[colData(example_sce)$cell_type=='b.cells',])[effect_gene_id]
# colMeans(example_para$mean_mat[colData(example_sce)$cell_type=='regulatory.t',])[effect_gene_id]

# colMeans(example_para$mean_mat[colData(example_sce)$cell_type=='b.cells',])[noneffect_gene_id]==colMeans(example_para$mean_mat[colData(example_sce)$cell_type=='regulatory.t',])[noneffect_gene_id]


# all_gene_index=seq(dim(example_sce)[1])-1
