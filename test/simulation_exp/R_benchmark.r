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
suppressPackageStartupMessages({
    library(optparse)
    library(PADOG, verbose = FALSE)
    library(SCPA, verbose = FALSE)
    library(SingleCellExperiment, verbose = FALSE)
    library(tidyverse, verbose = FALSE)
    library(reticulate, verbose = FALSE)
    library(zellkonverter, verbose = FALSE)
    library(scater, verbose = FALSE)
    library(msigdbr, verbose = FALSE)
    library(vctrs, verbose = FALSE)
    library(tibble, verbose = FALSE)
    library(limma, verbose = FALSE)
    library(optparse, verbose = FALSE)
    library(scuttle, verbose = FALSE)
    library(dplyr, verbose = FALSE)
    library(circlize, verbose = FALSE)
})




# 定义命令行参数
option_list <- list(
    make_option(
        c("--bench_param"),
        type = "character",
        default = NULL,
        help = "Benchmark axis to vary: overlap, sample, or sparsity , null_overlap, null_sparsity, beta",
        metavar = "character"
    ),
    make_option(
        c("--method"),
        type = "character",
        default = "all",
        help = "Which method to analyze (PADOG, CAMERA, SCPA, or all)",
        metavar = "character"
    ),
    make_option(
        c("--is_lognorm"),
        type = "logical",
        default = TRUE,
        help = "Whether to log-normalize counts before analysis",
        metavar = "logical"
    ),
    make_option(
        c("--fdr"),
        type = "double",
        default = 0.2,
        help = "FDR threshold [default= %default]",
        metavar = "number"
    )
)

# 解析参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


camera_test <- function(example_sce, fdr = 0.2) {
    set.seed(42)

    all_pathways <- metadata(example_sce)$all_pathways

    cell_type_matrix <- model.matrix(~ colData(example_sce)$cell_type)

    res <- camera(y = assays(example_sce)$logcounts, index = metadata(example_sce)$all_pathways, design = cell_type_matrix)

    selected_pathways <- rownames(res)[res$FDR < fdr]
    for (fdr in c(fdr)) {
        print(paste0("fdr=", fdr))
        print(paste0("CAMERA selects", length(rownames(res)[res$FDR < fdr]), "pathways"))

        real_power <- sum(metadata(example_sce)$effect_pathways %in% selected_pathways) / length(metadata(example_sce)$effect_pathways)
        real_fdr <- sum(!(selected_pathways %in% metadata(example_sce)$effect_pathways)) / length(selected_pathways)
        print(paste0("real_power=", real_power))
        print(paste0("real_fdr=", real_fdr))
    }
    return(c(real_fdr, real_power, length(selected_pathways)))
}


repeat_padog <- function(example_sce, fdr = 0.2) {
    cell_type <- recode(colData(example_sce)$cell_type, `0` = "c", `1` = "d")

    myr <- padog(
        esetm = logcounts(example_sce),
        group = cell_type,
        gslist = metadata(example_sce)$all_pathways,
        verbose = FALSE,
        Nmin = 3,
        NI = 50,
        dseed = 1,
        ncr = 100,
        parallel = TRUE
    )

    for (fdr in c(fdr)) {
        print(paste0("fdr=", fdr))
        print(paste0("select ", length(myr$Name[myr$Ppadog < fdr]), "pathways"))

        real_fdr <- 1 - sum(myr$Name[myr$Ppadog < fdr] %in% metadata(example_sce)$effect_pathways) / length(myr$Name[myr$Ppadog < fdr])
        print(paste0("real_fdr=", real_fdr))

        real_power <- sum(metadata(example_sce)$effect_pathways %in% myr$Name[myr$Ppadog < fdr]) / length(metadata(example_sce)$effect_pathways)
        print(paste0("real_power=", real_power))
    }

    return(c(real_fdr, real_power, length(myr$Name[myr$Ppadog < fdr])))
}


SCPA_test <- function(example_sce, fdr = 0.2) {
    pathways <- lapply(names(metadata(example_sce)$all_pathways), function(key) {
        tb <- data.frame(
            Pathway = rep(key, length(metadata(example_sce)$all_pathways[[key]])),
            Genes = unlist(metadata(example_sce)$all_pathways[[key]])
        )
        return(as_tibble(tb))
    })

    sce0 <- sce_extract(example_sce, assay_name = "logcounts", meta1 = "cell_type", value_meta1 = 0)
    sce1 <- sce_extract(example_sce, assay_name = "logcounts", meta1 = "cell_type", value_meta1 = 1)

    scpa_out <- compare_pathways(
        samples = list(sce0, sce1),
        pathways = pathways,
        parallel = TRUE,
        min_genes = 1
    )

    fdr_values <- c(fdr)
    for (fdr in fdr_values) {
        print(paste0("fdr=", fdr))
        selected_pathways <- scpa_out %>%
            filter(adjPval < fdr) %>%
            select(Pathway)

        real_power <- sum(metadata(example_sce)$effect_pathways %in% selected_pathways[, 1]) / length(metadata(example_sce)$effect_pathways)
        real_fdr <- sum(!(selected_pathways[, 1] %in% metadata(example_sce)$effect_pathways)) / nrow(selected_pathways)
        print(paste0("SCPA select ", length(selected_pathways[, 1]), " pathways"))
        print(paste0("real_power=", real_power))
        print(paste0("real_fdr=", real_fdr))

        print("-------------------------------------------------")
    }
    # print(metadata(example_sce)$effect_pathways)
    return(c(real_fdr, real_power, length(selected_pathways[, 1])))
}

seeds <- 46:75

if (opt$bench_param == "overlap") {
    params <- seq(0, 10, by = 2)
    param_prefix <- "n_overlap_"
    filename_builder <- function(path, param, seed) {
        paste0(path, "/simu_scRNAseq_100pathways_", param, "overlaplog_seed=", seed, ".h5ad")
    }
} else if (opt$bench_param == "sample") {
    params <- seq(600, 1800, by = 240)
    param_prefix <- "n_sample_"
    filename_builder <- function(path, param, seed) {
        paste0(path, "/simu_scRNAseq_100pathways_", param, "sample_seed=", seed, ".h5ad")
    }
} else if (opt$bench_param == "sparsity") {
    params <- seq(6, 16, by = 2)
    param_prefix <- "n_sparsity_"
    filename_builder <- function(path, param, seed) {
        paste0(path, "/simu_scRNAseq_100pathways_", param, "sparsity_seed=", seed, ".h5ad")
    }
} else if (opt$bench_param == "null_overlap") {
    params <- seq(0, 10, by = 2)
    param_prefix <- "n_overlap_"
    filename_builder <- function(path, param, seed) {
        paste0(path, "/simu_scRNAseq_100pathways_", param, "overlap_null_seed=", seed, ".h5ad")
    }
} else if (opt$bench_param == "null_sparsity") {
    params <- seq(6, 16, by = 2)
    param_prefix <- "n_sparsity_"
    filename_builder <- function(path, param, seed) {
        paste0(path, "/simu_scRNAseq_100pathways_", param, "sparsity_null_seed=", seed, ".h5ad")
    }
} else if (opt$bench_param == "beta") {
    params <- c(0, 0.5, 1, 1.5, 2, 3)
    param_prefix <- ""
    filename_builder <- function(path, param, seed) {
        paste0(path, "/simu_scRNAseq_100pathways_beta=", param, "_seed=", seed, ".h5ad")
    }
} else {
    stop("--bench_param must be one of overlap, sample, or sparsity")
}

fdr_list <- list()
power_list <- list()
n_selected_list <- list()


for (param in params) {
    print(paste0(opt$bench_param, "=", param))
    path <- file.path(
        "./data/simulation",
        paste0(opt$bench_param, "_data"),
        paste0(param_prefix, "beta", param)
    )
    if (opt$bench_param == "beta") {
        path <- paste0("./data/simulation/beta_data/beta=", param)
    }
    fdr_values <- numeric()
    power_values <- numeric()
    n_selected_values <- numeric()
    for (seed in seeds) {
        print(paste0("seed=", seed))

        file_name <- filename_builder(path, param, seed)
        example_sce <- readH5AD(file_name, reader = "R", use_hdf5 = FALSE)

        if (opt$is_lognorm) {
            example_sce <- logNormCounts(example_sce)
        }

        if (opt$method == "PADOG") {
            fdr_pow_nselected <- repeat_padog(example_sce, opt$fdr)
        } else if (opt$method == "SCPA") {
            fdr_pow_nselected <- SCPA_test(example_sce, opt$fdr)
        } else if (opt$method == "CAMERA") {
            fdr_pow_nselected <- camera_test(example_sce, opt$fdr)
        }

        fdr_values <- c(fdr_values, fdr_pow_nselected[1])
        power_values <- c(power_values, fdr_pow_nselected[2])
        n_selected_values <- c(n_selected_values, fdr_pow_nselected[3])
        print(fdr_pow_nselected)
        print("--------------------------------------------------------------------")
    }

    print(paste("FDR values:", fdr_values))
    print(paste("Power values:", power_values))
    print(paste("Number of selected pathways:", n_selected_values))

    fdr_list[[as.character(param)]] <- fdr_values
    power_list[[as.character(param)]] <- power_values
    n_selected_list[[as.character(param)]] <- n_selected_values
    print("===============================================")
}

fdr_df <- as.data.frame(fdr_list)
power_df <- as.data.frame(power_list)
n_selected_df <- as.data.frame(n_selected_list)

results_dir <- file.path(
    "./results/simulation_exp",
    paste0(opt$method, "_", opt$bench_param, "_results")
)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
    fdr_df,
    file.path(results_dir, paste0(opt$method, "_", opt$bench_param, "_fdr_values.csv")),
    row.names = FALSE
)

write.csv(
    power_df,
    file.path(results_dir, paste0(opt$method, "_", opt$bench_param, "_power_values.csv")),
    row.names = FALSE
)

write.csv(
    n_selected_df,
    file.path(results_dir, paste0(opt$method, "_", opt$bench_param, "_n_selected_values.csv")),
    row.names = FALSE
)
