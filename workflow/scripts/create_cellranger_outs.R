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
sample_id <- snakemake@wildcards[["sample"]]

all_df <- all_df %>% filter(sample == sample_id)

list_pools <- split(all_df$BARCODE, all_df$pool)

pools <- names(list_pools)

print(pools)

all_features <- NULL
for (pool in pools) {
    pools_bc <- list_pools[[pool]]

    pool_matrix_path <- file.path("results/assign_cells",pool,"outs/filtered_feature_bc_matrix/matrix.mtx.gz")
    sce <- read10xCounts(dirname(pool_matrix_path))
    all_features <- c(all_features,rowData(sce)$Symbol)
}
all_features <- unique(all_features)

all_sparse_matrix_filtered <- NULL
all_barcodes_filtered <- NULL

all_sparse_matrix <- NULL
all_barcodes <- NULL

for (pool in pools) {

    message("Processing pool ", pool)
    pools_bc <- list_pools[[pool]]

    pool_matrix_path <- file.path("results/assign_cells",pool,"outs/filtered_feature_bc_matrix/matrix.mtx.gz")

    #print(m_file)
    sce <- read10xCounts(dirname(pool_matrix_path))
    barcodes <- colData(sce)$Barcode
    new_barcodes <- paste0(colData(sce)$Barcode, "-", gsub("_","-",pool))
    missing_genes <- setdiff(all_features, rowData(sce)$Symbol)
    num_samples <- ncol(sce)
    # Create a new matrix filled with zeros for the missing genes
    new_genes_matrix <- matrix(
      0,
      nrow = length(missing_genes),
      ncol = num_samples,
      dimnames = list(missing_genes, colnames(sce))
    )
    sparse_matrix <- rbind(counts(sce), new_genes_matrix)
    
    sparse_matrix_filtered <- sparse_matrix[, which(barcodes %in% pools_bc)]
    barcodes_filtered <- new_barcodes[barcodes %in% pools_bc]

    all_sparse_matrix_filtered <- cbind(all_sparse_matrix_filtered, sparse_matrix_filtered)
    all_barcodes_filtered <- c(all_barcodes_filtered, barcodes_filtered)

    all_sparse_matrix <- cbind(all_sparse_matrix, sparse_matrix)
    all_barcodes <- c(all_barcodes, new_barcodes)
}

# Create the filtered directory
sample_dir_filtered <- paste0("results/cellranger_gex_out/", sample_id, "/outs/filtered_feature_bc_matrix/")
dir.create(sample_dir_filtered, recursive = T)
write10xCounts(path = sample_dir_filtered, 
                x = all_sparse_matrix_filtered, 
                barcodes = all_barcodes_filtered, 
                gene.id = all_features, 
                version = "3",
                overwrite = T)

# Create the raw directory
sample_dir_raw <- paste0("results/cellranger_gex_out/", sample_id, "/outs/raw_feature_bc_matrix/")
dir.create(sample_dir_raw, recursive = T)
write10xCounts(path = sample_dir_raw, 
                x = all_sparse_matrix, 
                barcodes = all_barcodes, 
                gene.id = all_features, 
                version = "3",
                overwrite = T)

file.create(snakemake@output[[1]])

sink()