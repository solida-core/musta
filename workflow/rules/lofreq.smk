rule lofreq:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/lofreq/{sample}_somatic_final_minus-dbsnp.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/lofreq/{sample}_somatic_final_minus-dbsnp.indels.vcf.gz",
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        dbsnp=config.get("resources").get("dbsnp"),
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/lofreq/{sample}_",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/rollcall/{sample}.lofreq.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/lofreq.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "lofreq somatic "
        "-f {params.genome} "
        "-t {input.tumoral} "
        "-n {input.normal} "
        "--threads {threads} "
        "-l {params.intervals} "
        "-d {params.dbsnp} "
        "-o {params.out} "
        ">& {log} "

