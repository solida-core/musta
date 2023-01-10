rule get_sample_names_vardict:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        normal=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/vardict/tmp/{sample}.normal.samplename.txt",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/vardict/tmp/{sample}.tumor.samplename.txt",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
    log:
        normal=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/vardict/{sample}.gsn.normal.log",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/vardict/{sample}.gsn.tumor.log",
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
        normal_name=rules.get_sample_names_vardict.output.normal,
        tumor_name=rules.get_sample_names_vardict.output.tumor,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/vardict/{sample}.somatic.vardict.vcf",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/vardict.yaml"
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        normal_name= lambda wildcards,input: get_name(input.normal_name),
        tumor_name=lambda wildcards, input: get_name(input.tumor_name),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "vardict "
        "-G {params.genome} "
        "-f 0.01 -N {params.tumor_name} "
        "-b '{input.tumoral} | {input.normal}' "
        "-c 1 -S 2 -E 3 "
        "{params.intervals} "
        "| testsomatic.R | var2vcf_paired.pl "
        "-N '{params.tumor_name}|{params.normal_name}' "
        "-f 0.01 > {output} "
