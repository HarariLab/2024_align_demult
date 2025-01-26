import os
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

rule filter_all:
    input: expand("results/filter_vcf/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule filter_vcf:
    input: "data/genotypes/{pool}.vcf.gz"
    params:
        MAF = config["MAF"],
        R2 = config["R2"]
    output: 
        done = "results/filter_vcf/{pool}.done",
        vcf = "results/filter_vcf/{pool}.filtered.vcf.gz"
    benchmark:
        "benchmarks/filter_vcf/{pool}.benchmark.txt"
    log: 
        "logs/filter_vcf/{pool}.out"
    conda: 
        config["conda_env"]
    shell:
        """
        bcftools filter -e 'INFO/MAF<{params.MAF}' {input} | bcftools filter \
        -e 'INFO/R2<{params.R2}' -Oz --output {output.vcf} && \
        touch {output.done}
        """

