import os
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

def get_vcfs(wildcards):
    bam = f"results/cellranger_arc_count/{wildcards.pool}/outs/atac_possorted_bam.bam"
    if(not(config.get("skip_vcf_sort", False))):
        vcf = f"results/updated_vcf/{wildcards.pool}.updated.sorted.vcf.gz"
        return { "bam": bam, "vcf": vcf }
    else:
        vcf = f"data/genotypes/{wildcards.pool}.vcf.gz"
        return { "bam": bam, "vcf": vcf }

rule  demuxlet_atac_all:
    input:
        expand("results/demuxlet_atac/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule demuxlet_atac:
    input:
        unpack(get_vcfs)
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
            --out results/demuxlet_atac/{wildcards.pool} \
            1> {log.out} \
            2> {log.err}
            
        touch {output}
        """
