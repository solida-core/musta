rule get_sample_names:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        normal=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "samplename/{sample}.normal.samplename.txt",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "samplename/{sample}.tumor.samplename.txt",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
    log:
        normal=resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "samplename/{sample}.gsn.normal.log",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "samplename/{sample}.gsn.tumor.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk "
        "--java-options {params.custom} "
        "GetSampleName "
        "-I {input.tumoral} "
        "-O {output.tumor} "
        ">& {log.tumor} && "
        "gatk "
        "--java-options {params.custom} "
        "GetSampleName "
        "-I {input.normal} "
        "-O {output.normal} "
        ">& {log.normal} "
