import os
configfile: "config/config.yaml"
include: "import_libraries.smk"

POOLS_LIST = POOLS.pool_id.unique().tolist()

rule sort_all:
    input: 
        expand("results/updated_vcf/{pool}.done", pool = POOLS_LIST)
    default_target: True

rule sort_vcf:
    input:
        bam = "results/cellranger_arc_count/{pool}/outs/gex_possorted_bam.bam",
        vcf = "data/genotypes/{pool}.vcf.gz"
    output: temp("results/sort_vcf/{pool}.sorted.vcf.gz")
    log: "logs/sort_vcf/{pool}.out"
    conda:
        config["conda_env"]
    shell:
        """
        picard SortVcf -INPUT {input.vcf} \
            -OUTPUT {output} \
            -R {input.bam} 2> {log}
        """

rule update_vcf_dict:
    input:
        bam = "results/cellranger_arc_count/{pool}/outs/gex_possorted_bam.bam",
        vcf = "results/sort_vcf/{pool}.sorted.vcf.gz"
    output: 
        done = "results/updated_vcf/{pool}.done",
        vcf = "results/updated_vcf/{pool}.updated.sorted.vcf.gz"
    log: 
        out = "logs/update_vcf_dict/{pool}.out",
        err = "logs/update_vcf_dict/{pool}.err"
    conda:
        config["conda_env"]
    shell:
        """
        picard UpdateVcfSequenceDictionary -INPUT {input.vcf} \
            -OUTPUT results/updated_vcf/{wildcards.pool}.updated.sorted.vcf.gz \
            -SD {input.bam} \
            1> {log.out} \
            2> {log.err}
        touch {output.done}
        """