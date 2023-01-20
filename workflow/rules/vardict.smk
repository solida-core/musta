rule vardict:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/vardict/{sample}.somatic.vardict.vcf.gz",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/vardict.yaml"
        ),
    params:
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/vardict/{sample}.somatic.vardict.vcf",
        ),
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
        "-f 0.01 > {params.out} ; "
        "bgzip -c {params.out} > {output} "
