rule MuSE_call:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),

    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/muse/{sample}.MuSE.txt"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/muse.yaml"
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),

        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/muse/{sample}",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/muse/{sample}.muse_call.log",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "MuSE call "
        "-f {params.genome} "
        "-O {params.out} "
        "{input.tumoral} " ## tumoral bam (positional)
        "{input.normal} " ## normal bam
        ">& {log} "

rule MuSE_sump:
    input:
        rules.MuSE_call.output,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
           "variant_calling/muse/{sample}.somatic.muse.vcf"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/muse.yaml"
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        dbsnp=config.get("resources").get("dbsnp"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/muse/{sample}.muse_sump.log",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "MuSE sump "
        "-I {input} " 
        "-D {params.dbsnp} "
        "-E "
        "-O {output} "
        ">& {log} "

rule MuSe_out:
    input:
        vcf=rules.MuSE_sump.output,
        normal_name= rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        resolve_results_filepath(
           config.get("paths").get("results_dir"),
            "variant_calling/muse/{sample}.somatic.muse.vcf.gz"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml"
        ),
    params:
        normal_name=lambda wildcards, input: get_name(input.normal_name),
        tumor_name=lambda wildcards, input: get_name(input.tumor_name),

    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "sed 's/NORMAL/{params.normal_name}/g' {input.vcf} ; "
        "sed 's/TUMOR/{params.tumor_name}/g' {input.vcf} ; "
        "bgzip -c {input.vcf} > {output} "
