# gen_sample.r with one change: DE genes are drawn from the own genes of a DE pathway (before the overlap step), see gen_overlap_own.r.
# Parameterised copy of simu_sample_diffsize.r for the revision (only seeds/values/cores changed).
# Usage: Rscript gen_sample.r <values comma-separated> <seed_from> <seed_to> <n_cores>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: Rscript generator.r values seed_from seed_to cores")
DATA_ROOT <- Sys.getenv("SCKP_SIM_DIR", "data/simulation")
REV_VALUES <- as.numeric(strsplit(args[1], ",")[[1]])
REV_SEEDS <- as.integer(args[2]):as.integer(args[3])
NCORES <- as.integer(args[4])
library(scDesign3, verbose = FALSE)
library(SingleCellExperiment)
library(DuoClustering2018)
library(tidyverse)
library(zellkonverter)
library(reticulate)
library(scuttle)
library(circlize)
library(GSA)
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
set.seed(123)



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

generate_sce_with_sample <- function(example_sce, n_sample, DE_pathways, groups, seed, save = TRUE) {
  path <- paste0(file.path(DATA_ROOT, "sampleown_data/sampleown="), n_sample)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  file_name <- paste0(path, "/simu_scRNAseq_100pathways_", n_sample, "sampleown_seed=", seed, ".h5ad")
  if (file.exists(file_name) && FALSE) {
    print(paste(file_name, "exists, skip simulation."))
    return(NULL)
  }

  # Add n_overlap genes to each pathway
  groups_own <- groups # gene sets before the overlap step
  groups_with_overlap <- groups
  n_overlap <- 5
  if (n_overlap != 0) {
    for (i in seq_along(groups)) {
      if (i != length(groups_with_overlap)) {
        # 获取下一个通路的基因
        next_pathway_genes <- groups_with_overlap[[i + 1]]
        n_borrow <- min(n_overlap, length(next_pathway_genes))

        previous_adjacent_genes <- groups_with_overlap[[i]][(length(groups_with_overlap[[i]]) - n_borrow + 1):length(groups_with_overlap[[i]])]

        adjacent_genes <- next_pathway_genes[1:n_borrow]

        if (n_borrow > 0) {
          groups_with_overlap[[i + 1]] <- unique(c(previous_adjacent_genes, groups_with_overlap[[i + 1]]))

          groups_with_overlap[[i]] <- unique(c(groups_with_overlap[[i]], adjacent_genes))
        }
      }
    }
  }

  groups <- groups_with_overlap
  # 给出两个通路之间的重叠基因示例
  if (n_overlap > 0) {
    for (i in 1:5) {
      overlap_genes <- intersect(groups[[i]], groups[[i + 1]])
      print(paste("Overlap between pathway", i, "and", i + 1, ":", paste(overlap_genes, collapse = ", ")))
    }
  }

  if (file.exists(file_name)) { print(paste("exists", file_name)); return(NULL) }
  set.seed(seed)
  example_data <- construct_data(
    sce = example_sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "1",
    ncell = n_sample,
  )

  example_marginal <- fit_marginal(
    data = example_data,
    predictor = "gene",
    mu_formula = "cell_type",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = NCORES,
    usebam = FALSE
  )

  example_copula <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = example_marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = NCORES,
    input_data = example_data$dat,
  )

  # --- 修改开始: 重构相关性矩阵 ---
  print("Modifying copula correlation matrix to enforce pathway structure...")

  # 1. 获取原始相关性矩阵 (假设 copula_list[[1]] 是矩阵，如果不是请调整访问方式)
  original_corr <- example_copula$copula_list[[1]]

  # 获取 copula 中实际存在的基因 (因为有部分基因可能被过滤掉了)
  valid_genes <- rownames(original_corr)
  # 2. 创建一个新的单位矩阵 (对角线为1，其余为0)
  new_corr <- diag(length(valid_genes))
  rownames(new_corr) <- valid_genes
  colnames(new_corr) <- valid_genes

  # 3. 遍历每个通路，保留通路内部基因的相关性
  # 注意：groups 中定义的基因可能有一部分不在valid_genes 里，需要取交集
  for (pathway_name in names(groups)) {
    # 获取该通路下的所有基因
    pathway_genes <- groups[[pathway_name]]

    # 只保留那些在 copula 矩阵中存在的基因
    genes_in_copula <- intersect(pathway_genes, valid_genes)

    # 如果该通路在过滤后至少剩下 2 个基因，则保留它们之间的相关性
    if (length(genes_in_copula) >= 2) {
      # 从原始矩阵中提取子块
      sub_matrix <- original_corr[genes_in_copula, genes_in_copula]

      # 填入新矩阵
      new_corr[genes_in_copula, genes_in_copula] <- sub_matrix
    }
  }

  # 4. 强制矩阵对称 (防止浮点误差) 并检查正定性
  new_corr <- (new_corr + t(new_corr)) / 2

  # [可选] 高斯 Copula 要求相关性矩阵必须是正定的。
  # 如果强行置零导致矩阵非正定，可能需要进行近似处理 (比如 nearPD)
  library(Matrix)
  new_corr <- as.matrix(Matrix::nearPD(new_corr, corr = TRUE)$mat)

  print(max(new_corr[upper.tri(new_corr)])) # 检查非对角线元素是否都被置零（或接近零）

  # 5. 将修改后的矩阵放回对象中
  example_copula$copula_list[[1]] <- new_corr
  print("Copula modification done.")

  example_para <- extract_para(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = example_marginal,
    n_cores = NCORES,
    family_use = "nb",
    new_covariate = example_data$newCovariate,
    data = example_data$dat
  )

  new_cell_type <- unlist(as.vector(example_data$newCovariate["cell_type"]))


  ground_truth <- rep(FALSE, length(groups)) # labels of pathways
  DEGs <- c() # list of DE genes
  for (pathway in DE_pathways) {
    # 选出 60% 的基因进行调控
    ground_truth[which(names(groups) == pathway)] <- TRUE
    effect_genes <- groups_own[pathway] # own genes only
    n <- length(unlist(effect_genes))
    sample_size <- ceiling(n * 0.6)
    effect_genes <- sample(unlist(effect_genes), sample_size)

    for (gene in effect_genes) {
      # 1. 随机决定是上调还是下调
      state <- sample(c(-1, 1), 1)

      # [改进点 1]: Log-Normal 分布生成 FC
      # 均值为 log(1.5) 或 log(2)，带有一定的波动 (sd=0.2)
      # 这样每个基因的 Fold Change 都是独特的，不再是均匀分布那么死板
      target_fc <- runif(1, 1.5, 2) # exp(0.5) ≈ 1.65倍

      # 找到两组细胞的索引
      idx_0 <- (new_cell_type == 0)
      idx_1 <- (new_cell_type == 1)

      # 获取基准均值（Control组的均值）
      # 注意：我们取 idx_0 的均值作为基准，保持 Control 组的稳定性
      base_vec <- example_para$mean_mat[idx_0, gene]
      base_mean <- mean(base_vec, na.rm = TRUE)

      # [改进点 2]: 引入细胞特异性噪声 (Cell-specific Jitter)
      # 我们不再让所有 y=1 的细胞均值都变成 base * fc
      # 而是让每个细胞的均值都在 target_fc 附近波动
      # 这种 "Jitter" 是打破基因间强共线性的关键！
      n_cells_1 <- sum(idx_1)
      cell_variation <- rnorm(n_cells_1, mean = 0, sd = 1) # 10% 的细胞间波动

      # 构建新的均值向量
      if (state == 1) {
        # 上调：Base * (FC + 细胞波动)
        # pmax 确保均值不为负
        new_means <- base_mean * (target_fc + cell_variation)
        example_para$mean_mat[idx_1, gene] <- pmax(0, new_means)
        # print(new_means)
      } else {
        # 下调：Base / (FC + 细胞波动)
        # 使用除法处理下调更符合生物学逻辑（1/FC）
        new_means <- base_mean / (target_fc + cell_variation)
        example_para$mean_mat[idx_1, gene] <- pmax(0, new_means)
        # print(new_means)
      }
      # Control 组保持不变 (或者你可以给 Control 组也加极其微小的扰动，如果需要的话)
      # example_para$mean_mat[idx_0, gene] <- ... (通常保持原样即可)
      DEGs <- c(DEGs, gene)
    }
  }


  # expression of non-DE pathways are set as the same
  non_DE_pathways <- setdiff(names(groups), DE_pathways)
  for (pathway in non_DE_pathways) {
    current_genes <- groups[[pathway]]
    # Ensure we don't accidentally revert a gene that was already modified (if overlap exists)
    current_genes <- setdiff(current_genes, DEGs)
    current_genes <- intersect(current_genes, colnames(example_para$mean_mat))

    n_cells_1 <- sum(idx_1)
    # cell_variation <- rnorm(n_cells_1, mean = 0, sd = 0.1) # 10%
    for (gene in current_genes) {
      idx_0 <- (new_cell_type == 0)
      idx_1 <- (new_cell_type == 1)
      base <- mean(example_para$mean_mat[idx_0, gene])
      example_para$mean_mat[idx_0, gene] <- base
      example_para$mean_mat[idx_1, gene] <- base
    }
  }

  effect_mat_values <- example_para$mean_mat[, DEGs[1]]
  print(paste0("Min: ", min(effect_mat_values), ", count: ", sum(effect_mat_values == min(effect_mat_values))))

  example_newcount <- simu_new(
    sce = example_sce,
    mean_mat = example_para$mean_mat,
    sigma_mat = example_para$sigma_mat,
    zero_mat = example_para$zero_mat,
    quantile_mat = NULL,
    copula_list = example_copula$copula_list,
    important_feature = NULL,
    n_cores = NCORES,
    family_use = "nb",
    input_data = example_data$dat,
    new_covariate = example_data$newCovariate,
    filtered_gene = example_data$filtered_gene
  )

  if (save) {
    new_sce <- SingleCellExperiment(
      assays = list(counts = example_newcount),
      colData = list(cell_type = as.vector(new_cell_type)),
      rowData = list(gene_name = rownames(example_sce)),
      metadata = list(effect_pathways = DE_pathways, all_pathways = groups, DEGs = DEGs)
    )
    print(dim(example_newcount))
    writeH5AD(new_sce, file_name)
    print(paste(file_name, "create success."))
  }
}

# --- 3. 生成固定规格的通路结构 (Fixed Size Distribution) ---
# 目标：40 个 15 基因, 20 个 20 基因, 20 个 25 基因
# 总基因需求：40*15 + 20*20 + 20*25 = 600 + 400 + 500 = 1500 基因
# 剩余基因将被删除

all_genes <- rowData(example_sce)$symbol

# 1. 定义通路大小向量
sizes_small <- rep(15, 40) # 40 个小通路
sizes_mid <- rep(20, 20) # 20 个中通路
sizes_large <- rep(25, 20) # 20 个大通路

# 合并并打乱顺序，或者保持有序以便观察
# 这里我们保持有序：前40个小，中间20个中，后20个大
pathway_sizes <- c(sizes_small, sizes_mid, sizes_large)
num_pathways <- length(pathway_sizes) # 应该 = 80
total_genes_needed <- sum(pathway_sizes) # 应该 = 1500

print(paste("Total pathways:", num_pathways))
print(paste("Total genes needed:", total_genes_needed))
print(paste("Original genes available:", length(all_genes)))

# 2. 筛选基因 (删除多余的基因)
if (length(all_genes) < total_genes_needed) {
  stop("Error: Not enough genes in example_sce to satisfy the requirement of 1500 genes.")
}

# 随机抽取 1500 个基因用于构建通路
genes_to_keep <- sample(all_genes, total_genes_needed)

# [Critical Step]: 更新 example_sce，只保留这 1500 个基因
# 这一步对应您要求的“删去剩下的 56 个基因”
example_sce <- example_sce[genes_to_keep, ]
print(paste("New example_sce dimensions:", paste(dim(example_sce), collapse = " x ")))

# 3. 构建 groups 列表
groups <- list()
start_idx <- 1
# 此时 example_sce 的行名已经只有 1500 个了，我们需要用 genes_to_keep 来切分
# 为了随机性，我们将 genes_to_keep 再打乱一次（虽然 sample 已经随机了，再做一次无妨）
# 或者直接按 sample 的顺序切分
current_genes_pool <- genes_to_keep
for (i in 1:num_pathways) {
  size <- pathway_sizes[i]
  end_idx <- start_idx + size - 1
  groups[[i]] <- current_genes_pool[start_idx:end_idx]
  start_idx <- end_idx + 1
}

# 4. 命名通路
# 清楚地标记大小，方便后续分析结果
# Pathway_Small_1...40, Pathway_Mid_1...20, Pathway_Large_1...20
names(groups) <- c(
  paste0("Pathway_Small_", 1:40),
  paste0("Pathway_Mid_", 1:20),
  paste0("Pathway_Large_", 1:20)
)
# --- 5. 运行模拟循环 ---

set.seed(42)

# [Config]: 选择 DE Pathways
# 策略：为了测试公平性，我们从三种规格中各选一些
# 选 4 个小的，3 个中的，3 个大的 = 10 个 DE Pathways
# Small (1-40), Mid (41-60), Large (61-80)

idx_small <- sample(1:40, 4)
idx_mid <- sample(41:60, 3)
idx_large <- sample(61:80, 3)

DE_pathways <- names(groups)[c(idx_small, idx_mid, idx_large)]

print("Selected DE Pathways:")
print(DE_pathways)
print("Sizes of DE Pathways:")
print(sapply(groups[DE_pathways], length))


# Extract regulatory T cells
selected_cells <- which(colData(example_sce)$phenoid == c("regulatory.t"))

example_sce <- example_sce[, selected_cells]
cell_id <- seq(1, dim(example_sce)[2])
sampled_cell_id <- sample(cell_id, size = length(cell_id) / 2, replace = FALSE)
cell_type <- rep(0, dim(example_sce)[2])
cell_type[sampled_cell_id] <- 1
colData(example_sce)$cell_type <- as.factor(cell_type)
example_sce

# Example usage:
# seeds <- 46:55
seeds <- REV_SEEDS
# n_samples <- c(600, 1200, 1800)
n_samples <- REV_VALUES
for (n_sample_value in n_samples) {
  for (seed in seeds) {
    generate_sce_with_sample(example_sce, n_sample_value, DE_pathways, groups, seed, save = TRUE)
  }
}
