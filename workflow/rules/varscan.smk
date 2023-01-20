rule pre_varscan2_tumoral:
    input:
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.tumoral.pileup",
        ),
    conda:
       resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "variant_calling/varscan/{sample}.tumoral.log",
        ),
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
        "-o {output} >& {log} "


rule pre_varscan2_germinal:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.normal.pileup",
        ),
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "variant_calling/varscan/{sample}.normal.log",
        ),
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
        "-o {output} >& {log}"


rule varscan2:
    input:
        tumoral=rules.pre_varscan2_tumoral.output,
        normal=rules.pre_varscan2_germinal.output,

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
            config.get("paths").get("log_dir"),
            "variant_calling/varscan/{sample}.varscan.log",
        ),
    shell:
        "varscan somatic "
        "{input.normal} " # normal samtools pileup
        "{input.tumoral} " # tumoral samtools pileup
        "--output-vcf "
        "--tumor-purity 0.2 "
        "--output-snp {output.snvs} "
        "--output-indel {output.indels} "
        ">& {log}  "


rule varscan2_out:
    input:
        snvs=rules.varscan2.output.snvs,
        indels=rules.varscan2.output.indels,
        normal_name= rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        snvs = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.somatic.varscan.snvs.vcf.gz",
        ),
        indels = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/varscan/{sample}.somatic.varscan.indels.vcf.gz",
        ),
    params:
        normal_name = lambda wildcards, input: get_name(input.normal_name),
        tumor_name = lambda wildcards, input: get_name(input.tumor_name),

    conda: resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "sed -i 's/NORMAL/{params.normal_name}/g' {input.snvs} ; "
        "sed -i 's/TUMOR/{params.tumor_name}/g' {input.snvs} ; "
        "sed -i 's/NORMAL/{params.normal_name}/g' {input.indels} ; "
        "sed -i 's/TUMOR/{params.tumor_name}/g' {input.indels} ; "
        "bgzip -c {input.snvs} > {output.snvs} ; "
        "bgzip -c {input.indels} > {output.indels} "
