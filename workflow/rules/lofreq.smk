rule lofreq:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        intervals=rules.prepare_bedfile.output.intervals,
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/lofreq/{sample}.somatic.lofreq.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/lofreq/{sample}.somatic.lofreq.indels.vcf.gz",
        ),
    params:
        genome=config.get("resources").get("reference"),
        dbsnp=config.get("resources").get("dbsnp"),
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/lofreq/{sample}.",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/lofreq/{sample}.lofreq.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/lofreq/{sample}.lofreq.txt",
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
        "-l {input.intervals} "
        "-d {params.dbsnp} "
        "-o {params.out} "
        "--call-indels "
        ">& {log} ; "
        "cp {params.out}somatic_final_minus-dbsnp.snvs.vcf.gz {output.snvs} && "
        "cp {params.out}somatic_final_minus-dbsnp.indels.vcf.gz {output.indels} "


rule lofreq_hold_on:
    input:
        snvs=rules.lofreq.output.snvs,
        indels=rules.lofreq.output.indels,
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.lofreq.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.lofreq.indels.vcf.gz",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "cp {input.snvs} {output.snvs} ; "
        "cp {input.indels} {output.indels} "
