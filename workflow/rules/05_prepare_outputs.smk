import os
import pandas as pd
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

rule create_directories_all:
    input: "results/cellranger_gex_out/create_cellranger_structure.done"
    default_target: True

rule map_all_singlets_and_doublets:
    input: expand("results/demuxlet_gex/{pool}.best", pool=POOLS_LIST)
    output: "results/cellranger_gex_out/all_singlets_and_doublets.csv"
    log: 
        err="logs/cellranger_gex_out/map_all_singlets_and_doublets.err",
        out="logs/cellranger_gex_out/map_all_singlets_and_doublets.out"
    conda: 
        config["conda_env"]
    script:
        "../scripts/prepare_all_cells_file.R"

# Function to read the CSV dynamically when needed
def get_samples_list(wildcards):
    df = pd.read_csv()
    return df["sample"].unique().tolist()

rule create_cellranger_structure:
    input: 
        all_df=rules.map_all_singlets_and_doublets.output,
        matrix=expand("results/assign_cells/{pool}/outs/filtered_feature_bc_matrix/matrix.mtx.gz", pool=POOLS_LIST),
    output:
        "results/cellranger_gex_out/create_cellranger_structure.done"
    log:
        err="logs/cellranger_gex_out/create_cellranger_structure.err",
        log="logs/cellranger_gex_out/create_cellranger_structure.log"
    conda:
        config["conda_env"]
    script:
        "../scripts/create_cellranger_outs.R"