configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

rule cellranger_arc_all:
    input:
        expand("results/cellranger_arc_count/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule create_libraries:
    input: config["libraries"]
    output: "results/create_libraries/{pool}.libraries.csv"
    run:
        # Create one library file by pool and keep them on the create_libraries diretory
        def create_libraries(libraries_path):
            for pool, sub_df in POOLS.groupby('pool_id'):
                sub_df.iloc[:,1:].to_csv(path_or_buf = f"results/create_libraries/{pool}.libraries.csv", index = False)
        create_libraries(input)
    
rule cellranger_arc_count:
    input:
        libraries = os.path.join("results/create_libraries", "{pool}.libraries.csv")
    params:
        cellranger_arc_path = config["cellrangerarc_path"],
        refdir = config["refdir"],
        jobmode = config["jobmode"]
    output: 
        "results/cellranger_arc_count/{pool}.done",
        raw_barcodes = "results/cellranger_arc_count/{pool}/outs/raw_feature_bc_matrix/barcodes.tsv.gz",
        barcodes = "results/cellranger_arc_count/{pool}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        bam_gex = "results/cellranger_arc_count/{pool}/outs/gex_possorted_bam.bam",
        bam_atac = "results/cellranger_arc_count/{pool}/outs/atac_possorted_bam.bam"
    benchmark:
        "benchmarks/cellranger_arc_count/{pool}.benchmark.txt"
    log:
       out = "logs/cellranger_arc_count/{pool}.out",
       err = "logs/cellranger_arc_count/{pool}.err"
    shell:
        """
        mkdir -p results/cellranger_arc_count
        cd results/cellranger_arc_count

        rm -rf {wildcards.pool} && \
        
        {params.cellranger_arc_path} count \
            --id={wildcards.pool} \
            --reference={params.refdir} \
            --libraries=../../{input.libraries} \
            --jobmode={params.jobmode} \
            --mempercore=3 \
            --maxjobs=12 \
            1> ../../{log.out} \
            2> ../../{log.err} && \
        cd ../.. && \
        touch {output}
        """