suppressPackageStartupMessages({
    library(tidyverse, verbose = FALSE)
    library(zellkonverter, verbose = FALSE)
    library(reticulate, verbose = FALSE)
    library(SingleCellExperiment, verbose = FALSE)
    library(scuttle)
    library(dplyr)
    library(GSA)
    library(optparse)
    library(PADOG)
})

# 定义命令行参数
option_list <- list(
    make_option(c("--data_file"),
        type = "character", default = NULL,
        help = "Path to the h5ad data file", metavar = "character"
    ),
    make_option(c("--cell_type"),
        type = "character", default = "all",
        help = "Which cell type to analyze", metavar = "character"
    ),
    make_option(c("--layer"),
        type = "character", default = NULL,
        help = "Data layer to use [default= %default]", metavar = "character"
    ),
    make_option(c("--gene_names"),
        type = "character", default = NULL, # 修复：None -> NULL
        help = "Gene names column [default= %default]", metavar = "character"
    ),
    make_option(c("--geneset_path"),
        type = "character", default = NULL,
        help = "Path to geneset file", metavar = "character"
    ),
    make_option(c("--y_cls"),
        type = "character", default = NULL,
        help = "Target column name", metavar = "character"
    ),
    make_option(c("--class1"),
        type = "character", default = NULL,
        help = "Control class name", metavar = "character"
    ),
    make_option(c("--class2"),
        type = "character", default = NULL,
        help = "Treatment class name", metavar = "character"
    ),
    make_option(c("--fdr"),
        type = "double", default = 0.2,
        help = "FDR threshold [default= %default]", metavar = "number"
    ),
    make_option(c("--method"),
        type = "character", default = "PADOG",
        help = "Analysis method: PADOG, CAMERA, or SCPA [default= %default]", metavar = "character"
    ),
    make_option(c("--gene_thresh"),
        type = "numeric", default = 15,
        help = "Number of least genes in a pathway", metavar = "number"
    ),
    make_option(c("--max_gene_thresh"),
        type = "numeric", default = 500,
        help = "Number of most genes in a pathway", metavar = "number"
    )
)

# 解析参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 修复：添加缺失的变量定义
method <- opt$method

# 参数验证
if (!method %in% c("PADOG", "CAMERA", "SCPA", "GSEA")) {
    stop("Method must be either 'PADOG' or 'CAMERA' or 'SCPA'")
}

# 必需参数检查
if (is.null(opt$data_file)) {
    stop("Data file must be specified")
}

if (is.null(opt$geneset_path)) {
    stop("Geneset path must be specified")
}

# 兼容读取 h5ad：优先使用 zellkonverter 的 R reader，失败时回退 python reader
safe_read_h5ad <- function(path) {
    message("[readH5AD] Trying reader='R' first for better compatibility...")
    sce <- tryCatch(
        {
            readH5AD(path, reader = "R")
        },
        error = function(e) {
            message("[readH5AD] reader='R' failed: ", conditionMessage(e))
            NULL
        }
    )

    if (!is.null(sce)) {
        message("[readH5AD] Successfully loaded with reader='R'.")
        return(sce)
    }

    message("[readH5AD] Falling back to reader='python'...")
    sce <- tryCatch(
        {
            readH5AD(path, reader = "python")
        },
        error = function(e) {
            msg <- paste0(
                "Failed to read h5ad with both readers.\n",
                "R-reader error has been printed above.\n",
                "Python-reader error: ", conditionMessage(e), "\n\n",
                "Likely cause: python anndata/h5py version incompatibility (e.g. IOSpec encoding_type='null').\n",
                "Recommended fixes:\n",
                "1) Update anndata in the environment used by zellkonverter/basilisk;\n",
                "2) Or pre-convert the h5ad with a newer scanpy/anndata and re-run."
            )
            stop(msg, call. = FALSE)
        }
    )

    message("[readH5AD] Successfully loaded with reader='python'.")
    sce
}

# methods
run_padog_analysis <- function(example_sce, layer, gene_name, y_cls, c_name, d_name, pathways, fdr, gene_thresh) {
    set.seed(123)

    if (!is.null(gene_name)) {
        if (!gene_name %in% colnames(rowData(example_sce))) {
            stop(paste0("[PADOG] gene_names column not found in rowData: ", gene_name), call. = FALSE)
        }
        rownames(example_sce) <- as.character(rowData(example_sce)[[gene_name]])
    }

    if (is.null(rownames(example_sce))) {
        stop("[PADOG] rownames(example_sce) is NULL. Please provide --gene_names with gene symbols.", call. = FALSE)
    }

    # Normalize and deduplicate feature names to maximize overlap with pathways.
    sce_gene_ids <- trimws(as.character(rownames(example_sce)))
    sce_gene_ids[is.na(sce_gene_ids) | sce_gene_ids == ""] <- paste0("NA_gene_", seq_len(sum(is.na(sce_gene_ids) | sce_gene_ids == "")))
    rownames(example_sce) <- make.unique(sce_gene_ids)

    if (!is.null(layer)) {
        X <- assays(example_sce)[[layer]]
    } else {
        X <- assays(example_sce)[[1]]
    }

    # padog 要求一个 base matrix，部分 h5ad 读取结果可能是 DelayedMatrix/dgCMatrix 结构
    if (!"matrix" %in% class(X)) {
        message("[PADOG] Converting assay to dense matrix (may increase memory usage)...")
        X <- as.matrix(X)
    }
    rownames(X) <- rownames(example_sce)
    print(dim(X))
    print(y_cls)

    clean_pathways <- lapply(pathways, function(genes) {
        gs <- trimws(as.character(genes))
        gs <- gs[!is.na(gs) & gs != ""]
        unique(gs)
    })

    filtered_pathways <- lapply(clean_pathways, function(genes) {
        intersect(genes, rownames(X))
    })
    filtered_pathways <- filtered_pathways[sapply(filtered_pathways, length) > gene_thresh & sapply(filtered_pathways, length) < opt$max_gene_thresh]

    overlap_n <- sum(rownames(X) %in% as.character(unlist(filtered_pathways)))
    message("[PADOG] Overlap genes before case harmonization: ", overlap_n)

    # Common issue: expression uses lower-case and GMT uses upper-case symbols.
    if (overlap_n <= 10) {
        x_upper <- toupper(rownames(X))
        filtered_pathways_upper <- lapply(clean_pathways, function(genes) {
            intersect(toupper(genes), x_upper)
        })
        filtered_pathways_upper <- filtered_pathways_upper[
            sapply(filtered_pathways_upper, length) > gene_thresh & sapply(filtered_pathways_upper, length) < opt$max_gene_thresh
        ]

        overlap_upper_n <- sum(x_upper %in% as.character(unlist(filtered_pathways_upper)))
        message("[PADOG] Overlap genes after case harmonization: ", overlap_upper_n)

        if (overlap_upper_n > overlap_n) {
            rownames(X) <- x_upper
            filtered_pathways <- filtered_pathways_upper
            overlap_n <- overlap_upper_n
        }
    }

    if (length(filtered_pathways) == 0) {
        stop("[PADOG] No pathways left after filtering. Try lowering --gene_thresh or checking gene ID type (symbol vs Ensembl).", call. = FALSE)
    }

    if (overlap_n <= 10) {
        stop(
            paste0(
                "[PADOG] Too few overlapping genes between expression and pathways (", overlap_n, "). ",
                "PADOG requires overlap > 10. ",
                "Please check that --gene_names uses gene symbols matching the GMT file, ",
                "or use a matching geneset database."
            ),
            call. = FALSE
        )
    }

    message("[PADOG] Filtered pathways: ", length(filtered_pathways))

    cell_type <- recode(colData(example_sce)[[y_cls]], !!c_name := "c", !!d_name := "d")
    cell_type <- factor(cell_type, levels = c("c", "d"))
    print(table(cell_type))

    myr <- padog(
        esetm = X,
        group = cell_type,
        gslist = filtered_pathways,
        verbose = FALSE,
        Nmin = gene_thresh,
        NI = 10,
        dseed = 1,
        ncr = 50,
        parallel = TRUE
    )

    # print and save results
    print(paste0("fdr=", fdr))
    print(paste0("select ", length(myr$Name[myr$Ppadog < fdr]), "pathways"))
    print(myr$Name[myr$Ppadog < fdr])

    return(myr$Name[myr$Ppadog < fdr])
}

run_camera_test <- function(example_sce, layer, gene_name, y_cls, c_name, d_name, pathways, fdr, gene_thresh) {
    set.seed(123)

    if (!is.null(gene_name)) {
        rownames(example_sce) <- rowData(example_sce)[[gene_name]]
    }

    if (!is.null(layer)) {
        X <- assays(example_sce)[[layer]]
    } else {
        X <- assays(example_sce)[[1]]
    }

    print(dim(X))

    # 修复：使用传入的pathways参数而不是未定义的all_pathways
    filtered_pathways <- lapply(pathways, function(genes) {
        intersect(genes, rownames(example_sce))
    })
    filtered_pathways <- filtered_pathways[sapply(filtered_pathways, length) > gene_thresh & sapply(filtered_pathways, length) < opt$max_gene_thresh]
    print(paste("Number of filtered pathways:", length(filtered_pathways)))

    print(y_cls)

    cell_type_matrix <- model.matrix(~ colData(example_sce)[[y_cls]])

    print(head(cell_type_matrix))

    res <- camera(y = X, index = filtered_pathways, design = cell_type_matrix)

    print(paste0("fdr=", fdr))
    print(paste0("CAMERA selects ", length(rownames(res)[res$FDR < fdr]), " pathways"))
    print(rownames(res)[res$FDR < fdr])

    return(rownames(res)[res$FDR < fdr])
}

run_SCPA_analysis <- function(example_sce, layer, gene_name, y_cls, c_name, d_name, pathways, fdr, gene_thresh) {
    set.seed(123)

    if (!is.null(gene_name)) {
        rownames(example_sce) <- rowData(example_sce)[[gene_name]]
    }

    if (is.null(layer)) {
        assays(example_sce)$lognorm <- assays(example_sce)[[1]]
    } else {
        assays(example_sce)$lognorm <- assays(example_sce)[[layer]]
    }

    sce0 <- sce_extract(example_sce,
        assay_name = "lognorm",
        meta1 = y_cls, value_meta1 = c_name # 修复：使用c_name而不是sce0_name
    )
    sce1 <- sce_extract(example_sce,
        assay_name = "lognorm",
        meta1 = y_cls, value_meta1 = d_name # 修复：使用d_name而不是sce1_name
    )
    scpa_out <- compare_pathways(
        samples = list(sce0, sce1),
        pathways = pathways, # 修复：使用pathways而不是gmt_files
        parallel = TRUE,
        cores = 30,
        min_genes = gene_thresh, # 修复：使用传入的gene_thresh
        max_genes = opt$max_gene_thresh,
        downsample = 1000
    )

    selected_pathways <- scpa_out %>%
        filter(adjPval < fdr) %>%
        select(Pathway)

    print(paste("SCPA selects", nrow(selected_pathways), "pathways")) # 修复：nrow而不是length

    return(selected_pathways)
}

run_GSEA_analysis <- function(example_sce, layer, gene_name, y_cls, c_name, d_name, pathways, fdr, gene_thresh) {
    set.seed(123)
    if (!is.null(gene_name)) {
        rownames(example_sce) <- rowData(example_sce)[[gene_name]]
    }

    if (!is.null(layer)) {
        X <- assays(example_sce)[[layer]]
    } else {
        X <- assays(example_sce)[[1]]
    }

    print(dim(X))

    # Filter pathways by minimum gene count
    filtered_pathways <- lapply(pathways, function(genes) {
        intersect(genes, rownames(example_sce))
    })

    filtered_pathways <- filtered_pathways[(sapply(filtered_pathways, length) > gene_thresh) & (sapply(filtered_pathways, length) < opt$max_gene_thresh)]
    print(paste("Number of filtered pathways:", length(filtered_pathways)))

    # Calculate log2 fold changes between classes
    class1_cells <- colData(example_sce)[[y_cls]] == c_name
    class2_cells <- colData(example_sce)[[y_cls]] == d_name

    # Calculate means and standard deviations for each class
    mean1 <- rowMeans(X[, class1_cells, drop = FALSE])
    mean2 <- rowMeans(X[, class2_cells, drop = FALSE])
    sd1 <- apply(X[, class1_cells, drop = FALSE], 1, sd)
    sd2 <- apply(X[, class2_cells, drop = FALSE], 1, sd)

    # Calculate signal-to-noise ratio
    gene_stats <- (mean2 - mean1) / (sd1 + sd2)
    names(gene_stats) <- rownames(X)

    # Remove NAs and infinite values
    gene_stats <- gene_stats[!is.na(gene_stats) & is.finite(gene_stats)]

    # Create ranked gene list
    gene_ranks <- sort(gene_stats, decreasing = TRUE)

    print(head(gene_ranks))

    # Run GSEA
    fgsea_results <- fgsea(
        pathways = filtered_pathways,
        stats = gene_ranks,
        minSize = gene_thresh,
        maxSize = 100,
        nperm = 1000
    )
    print(head(fgsea_results))

    # Filter significant pathways
    significant_pathways <- fgsea_results[fgsea_results$padj < fdr, ]$pathway

    print(paste0("fdr=", fdr))
    print(paste0("GSEA selects ", length(significant_pathways), " pathways"))
    print(significant_pathways)

    return(significant_pathways)
}


# Extract filename without extension and create corresponding folder
data_filename <- tools::file_path_sans_ext(basename(opt$data_file))
result_folder <- file.path("./results/real_exp", data_filename, paste0(opt$method, "_results"))
if (!dir.exists(result_folder)) {
    dir.create(result_folder, recursive = TRUE)
}

# read data and filter by cell type and disease type if specified
example_sce <- safe_read_h5ad(opt$data_file)
if (!is.null(opt$cell_type) && opt$cell_type != "all") {
    # 建议也加上 which() 以防 cell_type 列有 NA
    example_sce <- example_sce[, which(colData(example_sce)$cell_type == opt$cell_type)]
}
print(paste("After cell type filtering:", ncol(example_sce), "cells"))



# 修改后：使用 %in% 安全地处理 NA
target_classes <- c(opt$class1, opt$class2)
keep_cells <- colData(example_sce)[[opt$y_cls]] %in% target_classes
example_sce <- example_sce[, keep_cells]

colData(example_sce)[[opt$y_cls]] <- droplevels(colData(example_sce)[[opt$y_cls]])



# main process
if (method == "CAMERA") {
    library(limma)
    save_names <- paste0("CAMERA_", opt$class2, "_vs_", opt$class1, "_", opt$cell_type,"_", opt$y_cls, ".csv")
    save_names <- file.path(result_folder, save_names)

    # 读取pathway数据
    pathways <- GSA.read.gmt(opt$geneset_path)
    all_pathways <- setNames(pathways$genesets, pathways$geneset.names)

    selected_pathways <- run_camera_test(example_sce, opt$layer, opt$gene_names, opt$y_cls, opt$class1, opt$class2, all_pathways, opt$fdr, opt$gene_thresh)

    write.csv(data.frame(Pathway = selected_pathways), save_names, row.names = FALSE)
    print("save success")
}

if (method == "PADOG") {
    library(PADOG)
    save_names <- paste0("PADOG_", opt$class2, "_vs_", opt$class1, "_", opt$cell_type,"_", opt$y_cls, ".csv")
    save_names <- file.path(result_folder, save_names)

    # 读取pathway数据
    pathways <- purrr::quietly(GSA.read.gmt)(opt$geneset_path)$result
    all_pathways <- setNames(pathways$genesets, pathways$geneset.names)

    # 去除大于 100 genes 的 pathway


    selected_pathways <- run_padog_analysis(example_sce, opt$layer, opt$gene_names, opt$y_cls, opt$class1, opt$class2, all_pathways, opt$fdr, opt$gene_thresh)

    write.csv(data.frame(Pathway = selected_pathways), save_names, row.names = FALSE)
    print("save success")
}

if (method == "SCPA") {
    library(SCPA)
    save_names <- paste0("SCPA_", opt$class2, "_vs_", opt$class1, "_", opt$cell_type,"_", opt$y_cls, ".csv")
    save_names <- file.path(result_folder, save_names)

    selected_pathways <- run_SCPA_analysis(example_sce, opt$layer, opt$gene_names, opt$y_cls, opt$class1, opt$class2, opt$geneset_path, opt$fdr, opt$gene_thresh)

    write.csv(selected_pathways, save_names, row.names = FALSE)
    print("save success")
}

if (method == "GSEA") {
    library(fgsea)
    save_names <- paste0("GSEA_", opt$class2, "_vs_", opt$class1, "_", opt$cell_type,"_", opt$y_cls, ".csv")
    save_names <- file.path(result_folder, save_names)

    # 读取pathway数据
    pathways <- GSA.read.gmt(opt$geneset_path)
    all_pathways <- setNames(pathways$genesets, pathways$geneset.names)

    selected_pathways <- run_GSEA_analysis(example_sce, opt$layer, opt$gene_names, opt$y_cls, opt$class1, opt$class2, all_pathways, opt$fdr, opt$gene_thresh)

    write.csv(data.frame(Pathway = selected_pathways), save_names, row.names = FALSE)
    print("save success")
}
