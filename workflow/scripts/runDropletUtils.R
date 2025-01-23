## Capture the output and error logs
output_log <- file(snakemake@log[["out"]], open="wt")
error_log <- file(snakemake@log[["err"]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages({
    library(DropletUtils)
    library(dplyr)
})
set.seed(123)

RAW_COUNTS_DIR <- dirname(snakemake@input[[1]])
OUTPUT_DIR <- dirname(snakemake@output[["matrix"]])
sce <- read10xCounts(SAMPLE_DIR)
colnames(sce) <- colData(sce)$Barcode
sce <- sce[grepl("^ENSG",rownames(sce)),]
e.out <- emptyDrops(counts(sce))
keep_cells <- rownames(e.out %>% as.data.frame() %>% filter(!is.na(FDR) & e.out$FDR <= 0.01))
sce_filtered <- sce[,colnames(sce) %in% keep_cells]
dir.create(OUT_DIR, showWarnings = F, recursive = T)
write10xCounts(x = counts(sce_filtered),
            overwrite = TRUE,
            version = "3",
            path = OUT_DIR)
file.create(snakemake@output[[1]])


sink()
sink()


