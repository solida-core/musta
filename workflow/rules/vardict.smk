rule get_sample_names:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        normal=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/tmp/{sample}_normal.samplename.txt",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/tmp/{sample}_tumor.samplename.txt",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
    log:
        normal=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/gatk/getsamplename/{sample}.gsn.log",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/gatk/getsamplename/{sample}.gsn_tumor.log",
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


rule vardict:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/tmp/{sample}_normal.samplename.txt",
        ),
        tumor_name=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/tmp/{sample}_tumor.samplename.txt",
        ),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/vardict/{sample}.vardict.vcf",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/vardict.yaml"
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),

    shell:
        "vardict "
        "-G {params.genome} "
        "-f 0.01 -N {input.tumor_name} "
        "-b '{input.tumoral} | {input.normal}' "
        "-c 1 -S 2 -E 3 "
        "{params.intervals} "
        "| testsomatic.R | var2vcf_paired.pl "
        "-N '{input.tumor_name}|{input.normal_name}' "
        "-f 0.01 > {output} "
