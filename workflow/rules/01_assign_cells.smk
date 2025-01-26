import os
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

rule assign_cells_all:
    input: expand("results/assign_cells/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule run_droplet_utils:
    input:
        "results/cellranger_arc_count/{pool}/outs/raw_feature_bc_matrix/barcodes.tsv.gz"
    conda: 
        config["conda_env"]
    log:
        out = "logs/assign_cells/{pool}.out",
        err = "logs/assign_cells/{pool}.err"
    output:
        "results/assign_cells/{pool}.done",
        matrix = "results/assign_cells/{pool}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        features = "results/assign_cells/{pool}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        barcodes = "results/assign_cells/{pool}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    benchmark:
        "benchmarks/assign_cells/{pool}.benchmark.txt"
    script:
        "../scripts/runDropletUtils.R"
