import os
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

rule  demuxlet_atac_all:
    input:
        expand("results/demuxlet_atac/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule demuxlet_atac:
    input:
        bam = "results/cellranger_arc_count/{pool}/outs/atac_possorted_bam.bam",
        vcf = "results/updated_vcf/{pool}.updated.sorted.vcf.gz",
        barcodes = "results/cellranger_arc_count/{pool}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output: "results/demuxlet_atac/{pool}.done"
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
            --field GT --alpha 0 --alpha 0.5 --min-MQ 255 \
            --group-list {input.barcodes} \
            --out results/demuxlet_atac/{wildcards.pool} \
            1> {log.out} \
            2> {log.err}
            
        touch {output}
        """
