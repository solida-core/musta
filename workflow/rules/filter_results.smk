
rule filter_mutect_tumoronly:
    input:
        vcf=rules.mutect_tumoronly.output.vcf,
        cont_tab=rules.calculate_contamination_tumoronly.output.table,
        orient=rules.orientation_model_tumoronly.output,
        segments=rules.calculate_contamination_tumoronly.output.segment,
        stats=rules.mutect_tumoronly.output.stats,
    output:
        vcf=resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/tumoronly/{sample}.somatic.filtered.vcf.gz",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/detection/mutect/{sample}.filter_info.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk FilterMutectCalls "
        "--java-options {params.custom} "
        "-V {input.vcf} "
        "--tumor-segmentation {input.segments} "
        "--contamination-table {input.cont_tab} "
        "--ob-priors {input.orient} "
        "-O {output.vcf} "
        "--stats {input.stats} "
        "-R {params.genome} "
        "-L {params.intervals} "
        ">& {log} "


