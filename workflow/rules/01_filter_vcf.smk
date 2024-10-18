import os
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

rule filter_all:
    input: expand("results/filter_vcf/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule filter_vcf:
    input: "data/genotypes/{pool}.vcf.gz"
    output: 
        done = "results/filter_vcf/{pool}.done",
        vcf = "results/filter_vcf/{pool}.filtered.vcf.gz"
    log: "logs/filter_vcf/{pool}.out"
    conda: 
        config["conda_env"]
    shell:
        """
        bcftools filter --include 'MAF>=0.05' -Oz --output {output.vcf} {input}
        touch {output.done}
        """

