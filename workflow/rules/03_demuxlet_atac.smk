import os
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

def get_vcfs(wildcards):
    bam = f"results/cellranger_arc_count/{wildcards.pool}/outs/atac_possorted_bam.bam"
    barcodes = f"results/cellranger_arc_count/{wildcards.pool}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    if(not(config.get("skip_vcf_sort", False))):
        vcf = f"results/updated_vcf/{wildcards.pool}.updated.sorted.vcf.gz"
        return { "bam": bam, "vcf": vcf , "barcodes": barcodes }
    else:
        vcf = f"data/genotypes/{wildcards.pool}.vcf.gz"
        return { "bam": bam, "vcf": vcf, "barcodes": barcodes}

rule  demuxlet_atac_all:
    input:
        expand("results/demuxlet_atac/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule demuxlet_atac:
    input:
        unpack(get_vcfs)
    output: 
        "results/demuxlet_atac/{pool}.done"
    params: 
        cap_bq = config["cap_bq"],
        min_bq = config["min_bq"],
        min_mq = config["min_mq"],
        min_td = config["min_td"],
        excl_flag = config["excl_flag"],
        geno_error = config["geno_error"]
    log:
        out = "logs/demuxlet_atac/{pool}.out",
        err = "logs/demuxlet_atac/{pool}.err"
    benchmark:
        "benchmarks/demuxlet_atac/{pool}.benchmark.txt"
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
            1> {log.out} \
            2> {log.err}

        touch {output}
        """
