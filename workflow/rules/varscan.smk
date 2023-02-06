rule pre_varscan2_tumoral:
    input:
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/varscan/{sample}.tumoral.pileup",
        ),
    conda:
       resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/varscan/{sample}.tumoral.log",
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
            "detection/varscan/{sample}.normal.pileup",
        ),
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/varscan/{sample}.normal.log",
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
            "detection/varscan/{sample}.somatic.varscan.snvs.vcf",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/varscan/{sample}.somatic.varscan.indels.vcf",
        )
    conda: resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/varscan.yaml")
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/varscan/{sample}.varscan.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/varscan/{sample}.varscan.txt",
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
            "detection/varscan/{sample}.somatic.varscan.snvs.vcf.gz",
        ),
        indels = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/varscan/{sample}.somatic.varscan.indels.vcf.gz",
        ),
    params:
        normal_name = lambda wildcards, input: get_name(input.normal_name),
        tumor_name = lambda wildcards, input: get_name(input.tumor_name),

    conda: resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/tabix.yaml"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "sed -i 's/NORMAL/{params.normal_name}/g' {input.snvs} ; "
        "sed -i 's/TUMOR/{params.tumor_name}/g' {input.snvs} ; "
        "sed -i 's/NORMAL/{params.normal_name}/g' {input.indels} ; "
        "sed -i 's/TUMOR/{params.tumor_name}/g' {input.indels} ; "
        "bgzip -c {input.snvs} > {output.snvs} ; "
        "bgzip -c {input.indels} > {output.indels} "

rule varscan_hold_on:
    input:
        snvs=rules.varscan2_out.output.snvs,
        indels=rules.varscan2_out.output.indels,

    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.varscan.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.varscan.indels.vcf.gz",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "cp {input.snvs} {output.snvs} ; "
        "cp {input.indels} {output.indels} "

rule varscan_tbi:
    input:
        snvs=rules.varscan_hold_on.output.snvs,
        indels=rules.varscan_hold_on.output.indels,
    output:
        snvs = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.varscan.snvs.vcf.gz.tbi",
        ),
        indels = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.varscan.indels.vcf.gz.tbi",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/tabix.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "tabix -p vcf {input.snvs} ; "
        "tabix -p vcf {input.indels}  "
