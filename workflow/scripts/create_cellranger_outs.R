## Capture the output and error logs
output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages({
    library(DropletUtils)
    library(dplyr)
    library(Seurat)
})

all_df <- read.csv(snakemake@input[["all_df"]])
list_samples <- split(all_df$BARCODE, all_df$sample)

df_samples <- all_df %>% select(sample, pool) %>% unique()
pools <- setNames(df_samples$pool, df_samples$sample)

samples <- names(list_samples)

all_matrix_input <- snakemake@input[["matrix"]]
print(all_matrix_input)

for (s in samples) {

    message("Processing sample ", s)
    samples_bc <- list_samples[[s]]

    m_file <- all_matrix_input[grep(pools[s], all_matrix_input)]
    #print(m_file)
    sparse_matrix <- read10xCounts(dirname(m_file))
    sparse_matrix <- counts(sparse_matrix)
    sparse_matrix <- sparse_matrix[, which(colnames(sparse_matrix) %in% samples_bc)]

    features <- rownames(sparse_matrix)
    barcodes <- colnames(sparse_matrix)

    sample_dir <- paste0("results/cellranger_gex_out/", s, "/outs/filtered_feature_bc_matrix/")
    dir.create(sample_dir, recursive = T)
    write10xCounts(path = sample_dir, 
                   x = sparse_matrix, barcodes = barcodes, gene.id = features, 
                   version = "3",
                   overwrite = T)

}

file.create(snakemake@output[[1]])

sink()
sink()