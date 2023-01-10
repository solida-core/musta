rule get_sample_names_varscan:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        normal=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/tmp/{sample}.normal.samplename.txt",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/tmp/{sample}.tumor.samplename.txt",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
    log:
        normal=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/varscan/{sample}.gsn.normal.log",
        ),
        tumor=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/varscan/{sample}.gsn.tumor.log",
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


rule pre_varscan2_tumoral:
    input:
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=rules.get_sample_names_varscan.output.normal,
        tumor_name=rules.get_sample_names_varscan.output.tumor,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.tumoral.pileup",
        ),
    conda:
       resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),

    shell:
        "samtools "
        "mpileup "
        "-f {params.genome} "
        "-x -C 50 -q 40 -Q 30 "
        "-l {params.intervals} "
        "{input.tumoral} "
        "-o {output} "


rule pre_varscan2_germinal:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.normal.pileup",
        ),
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/samtools.yaml")

    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "samtools "
        "mpileup "
        "-f {params.genome} "
        "-x -C 50 -q 40 -Q 30 "
        "-l {params.intervals} "
        "{input.normal} "
        "-o {output} "


rule varscan2:
    input:
        tumoral=rules.pre_varscan2_tumoral.output,
        normal=rules.pre_varscan2_germinal.output
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.somatic.varscan.snvs.vcf",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.somatic.varscan.indels.vcf",
        )
    conda: resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/varscan.yaml")
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/varscan/{sample}.varscan.log",
        ),
    shell:
        "varscan somatic "
        "{input.normal} " # normal samtools pileup
        "{input.tumoral} " # tumoral samtools pileup
        "--output-vcf "
        "--tumor-purity 0.2 "
        "--output-snp {output.snvs} "
        "--output-indel {output.indels}"
        ">& {log}"
