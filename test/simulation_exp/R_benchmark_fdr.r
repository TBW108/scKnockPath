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
    )
)

# 解析参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

camera_test <- function(example_sce, fdr_values) {
    set.seed(42)
    cell_type_matrix <- model.matrix(~ colData(example_sce)$cell_type)
    res <- camera(
        y = assays(example_sce)$logcounts,
        index = metadata(example_sce)$all_pathways,
        design = cell_type_matrix
    )

    effect_pathways <- metadata(example_sce)$effect_pathways
    real_fdrs <- numeric(length(fdr_values))
    real_powers <- numeric(length(fdr_values))

    for (i in seq_along(fdr_values)) {
        fdr <- fdr_values[i]
        selected_pathways <- rownames(res)[res$FDR < fdr]

        print(paste0("fdr=", fdr))
        print(paste0("CAMERA selects", length(selected_pathways), "pathways"))

        real_powers[i] <- sum(effect_pathways %in% selected_pathways) / length(effect_pathways)
        real_fdrs[i] <- if (length(selected_pathways) == 0) {
            0
        } else {
            sum(!(selected_pathways %in% effect_pathways)) / length(selected_pathways)
        }

        print(paste0("real_power=", real_powers[i]))
        print(paste0("real_fdr=", real_fdrs[i]))
    }

    tibble(
        threshold = fdr_values,
        real_fdr = real_fdrs,
        real_power = real_powers
    )
}


repeat_padog <- function(example_sce, fdr_values) {
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

    effect_pathways <- metadata(example_sce)$effect_pathways
    real_fdrs <- numeric(length(fdr_values))
    real_powers <- numeric(length(fdr_values))

    for (i in seq_along(fdr_values)) {
        fdr <- fdr_values[i]
        selected <- myr$Name[myr$Ppadog < fdr]

        print(paste0("fdr=", fdr))
        print(paste0("select ", length(selected), " pathways"))

        real_fdrs[i] <- if (length(selected) == 0) {
            0
        } else {
            1 - sum(selected %in% effect_pathways) / length(selected)
        }
        real_powers[i] <- sum(effect_pathways %in% selected) / length(effect_pathways)

        print(paste0("real_fdr=", real_fdrs[i]))
        print(paste0("real_power=", real_powers[i]))
    }

    tibble(
        threshold = fdr_values,
        real_fdr = real_fdrs,
        real_power = real_powers
    )
}


SCPA_test <- function(example_sce, fdr_values) {
    pathways <- lapply(names(metadata(example_sce)$all_pathways), function(key) {
        tb <- data.frame(
            Pathway = rep(key, length(metadata(example_sce)$all_pathways[[key]])),
            Genes = unlist(metadata(example_sce)$all_pathways[[key]])
        )
        as_tibble(tb)
    })

    sce0 <- sce_extract(example_sce, assay_name = "logcounts", meta1 = "cell_type", value_meta1 = 0)
    sce1 <- sce_extract(example_sce, assay_name = "logcounts", meta1 = "cell_type", value_meta1 = 1)

    scpa_out <- compare_pathways(
        samples = list(sce0, sce1),
        pathways = pathways,
        parallel = TRUE,
        downsample = 1500,
        min_genes = 1
    )

    effect_pathways <- metadata(example_sce)$effect_pathways
    real_fdrs <- numeric(length(fdr_values))
    real_powers <- numeric(length(fdr_values))

    for (i in seq_along(fdr_values)) {
        fdr <- fdr_values[i]
        print(paste0("fdr=", fdr))

        selected_pathways <- scpa_out %>%
            filter(adjPval < fdr) %>%
            pull(Pathway)

        real_powers[i] <- sum(effect_pathways %in% selected_pathways) / length(effect_pathways)
        real_fdrs[i] <- if (length(selected_pathways) == 0) {
            0
        } else {
            sum(!(selected_pathways %in% effect_pathways)) / length(selected_pathways)
        }

        print(paste0("SCPA select ", length(selected_pathways), " pathways"))
        print(paste0("real_power=", real_powers[i]))
        print(paste0("real_fdr=", real_fdrs[i]))
        print("-------------------------------------------------")
    }

    tibble(
        threshold = fdr_values,
        real_fdr = real_fdrs,
        real_power = real_powers
    )
}

seeds <- 46:75
fdrs <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3)
fdr_keys <- as.character(fdrs)
fdr_list <- setNames(vector("list", length(fdrs)), fdr_keys)
power_list <- setNames(vector("list", length(fdrs)), fdr_keys)
n_overlap <- 4

for (seed in seeds) {
    print(paste0("seed=", seed))

    file_name <- paste0(
        "/disk/tanbowen/scKnockPath/data/simulation/overlap_data/n_overlap_beta",
        n_overlap,
        "/simu_scRNAseq_100pathways_", n_overlap,
        "overlaplog_seed=", seed, ".h5ad"
    )
    example_sce <- readH5AD(file_name, reader = "R", use_hdf5 = FALSE)

    if (opt$is_lognorm) {
        example_sce <- logNormCounts(example_sce)
    }

    metrics <- if (opt$method == "PADOG") {
        repeat_padog(example_sce, fdrs)
    } else if (opt$method == "SCPA") {
        SCPA_test(example_sce, fdrs)
    } else if (opt$method == "CAMERA") {
        camera_test(example_sce, fdrs)
    } else {
        stop("Unsupported method.")
    }

    print(metrics)

    for (i in seq_along(fdrs)) {
        key <- fdr_keys[i]
        fdr_list[[key]] <- c(fdr_list[[key]], metrics$real_fdr[i])
        power_list[[key]] <- c(power_list[[key]], metrics$real_power[i])
    }

    print("--------------------------------------------------------------------")
}

for (key in fdr_keys) {
    print(paste0("target_fdr=", key))
    print(fdr_list[[key]])
    print(power_list[[key]])
    print("===============================================")
}

fdr_df <- as.data.frame(fdr_list)
power_df <- as.data.frame(power_list)

results_dir <- file.path(
    "./results/simulation_exp",
    paste0(opt$method, "_fdrs_results")
)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
    fdr_df,
    file.path(results_dir, paste0(opt$method, "_fdrs_fdr_values.csv")),
    row.names = FALSE
)

write.csv(
    power_df,
    file.path(results_dir, paste0(opt$method, "_fdrs_power_values.csv")),
    row.names = FALSE
)
