configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

def get_vcfs(wildcards):
    bam = f"results/cellranger_arc_count/{wildcards.pool}/outs/gex_possorted_bam.bam"
    barcodes = f"results/assign_cells/{wildcards.pool}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    vcf = ""
    if(not(config.get("skip_vcf_filter", False))):
        vcf = f"results/filter_vcf/{wildcards.pool}.filtered.vcf.gz"
    else:
        vcf = f"data/genotypes/{wildcards.pool}.vcf.gz"
    return { "bam": bam, "vcf": vcf, "barcodes": barcodes }

rule demuxlet_gex_all:
    input:
        expand("results/demuxlet_gex/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule demuxlet_gex:
    input:
        unpack(get_vcfs)
    output: 
        "results/demuxlet_gex/{pool}.done"
    params: 
        cap_bq = config["cap_bq"],
        min_bq = config["min_bq"],
        min_mq = config["min_mq"],
        min_td = config["min_td"],
        excl_flag = config["excl_flag"],
        geno_error = config["geno_error"]
    log: "logs/demuxlet_gex/{pool}.out"
    benchmark:
        "benchmarks/demuxlet_gex/{pool}.benchmark.txt"
    conda:
        config["conda_env"]       
    shell:
        """
        demuxlet --sam {input.bam} --vcf {input.vcf} \
            --field GT --alpha 0 --alpha 0.5 \
            --geno-error {params.geno_error} \
            --min-BQ {params.min_bq} --cap-BQ {params.cap_bq} \
            --min-MQ {params.min_mq} --min-TD {params.min_td} \
            --excl-flag {params.excl_flag} \
            --group-list {input.barcodes} \
            --out results/demuxlet_gex/{wildcards.pool} \
            > {log} 2>&1 && \
        touch {output}
        """
