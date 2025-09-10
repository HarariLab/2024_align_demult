## Capture the output and error logs
output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages({
    library(DropletUtils)
    library(dplyr)
})

all_df <- read.csv(snakemake@input[["all_df"]])
pool_matrix_path <- snakemake@input[["matrix"]]

pool_path <- gsub(file.path("","outs","filtered_feature_bc_matrix","matrix.mtx.gz"),"",pool_matrix_path)
pool_id <- tail(unlist(strsplit(pool_path, split = "\\/")), n = 1)

all_df <- all_df %>% filter(pool == pool_id)

list_samples <- split(all_df$BARCODE, all_df$sample)

samples <- names(list_samples)

print(samples)

for (s in samples) {

    message("Processing sample ", s)
    samples_bc <- list_samples[[s]]

    m_file <- pool_matrix_path
    #print(m_file)
    sce <- read10xCounts(dirname(m_file))
    sce_filtered <- sce[, which(colData(sce)$Barcode %in% samples_bc) ]

    sparse_matrix <- counts(sce)
    sparse_matrix_filtered <- counts(sce_filtered)

    features <- rowData(sce)$Symbol
    barcodes <- colData(sce)$Barcode

    features_filtered <- rowData(sce_filtered)$Symbol
    barcodes_filtered <- colData(sce_filtered)$Barcode

    # Create the filtered directory
    sample_dir_filtered <- paste0("results/cellranger_gex_out/", s, "/outs/filtered_feature_bc_matrix/")
    dir.create(sample_dir_filtered, recursive = T)
    write10xCounts(path = sample_dir_filtered, 
                   x = sparse_matrix_filtered, barcodes = barcodes_filtered, gene.id = features_filtered, 
                   version = "3",
                   overwrite = T)

    # Create the raw directory
    sample_dir_raw <- paste0("results/cellranger_gex_out/", s, "/outs/raw_feature_bc_matrix/")
    dir.create(sample_dir_raw, recursive = T)
    write10xCounts(path = sample_dir_raw, 
                   x = sparse_matrix, barcodes = barcodes, gene.id = features, 
                   version = "3",
                   overwrite = T)

    

}

file.create(snakemake@output[[1]])

sink()
sink()