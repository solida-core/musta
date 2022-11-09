rule lofreq:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/call/{sample}_somatic_final_minus-dbsnp.snvs.vcf.gz",
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=resolve_single_filepath(
            config.get("paths").get("workdir"),
            resolve_single_filepath("resources", config.get("resources").get("bed")),
        ),
        param=config.get("params").get("gatk").get("Mutect"),
        dbsnp=config.get("params").get("lofreq").get("dbsnp"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/call/{sample}.lofreq.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/lofreq.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "lofreq "
        "somatic "
        "-f {params.genome} "
        "-t {input.tumoral} "
        "-n {input.normal} "
        "-l {params.intervals} "
        "-d {params.dbsnp} "
        "-o {sample}_ "
        ">& {log} "

