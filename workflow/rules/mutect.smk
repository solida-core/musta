rule mutect_matched:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        vcf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.somatic.vcf.gz",
        ),
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.tumor_normal.bam",
        ),
        fir=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.tumor_normal.f1r2.tar.gz",
        ),
        stats=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.somatic.vcf.gz.stats",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        param=config.get("params").get("gatk").get("Mutect"),
        germline_resource=config.get("params").get("gatk").get("germline"),
        normal_bam=lambda wildcards, input: get_name(input.normal_name),
        tumor_bam=lambda wildcards, input: get_name(input.tumor_name),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/{sample}.mutect.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/mutect/{sample}.mutect.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk "
        "--java-options {params.custom} "
        "Mutect2 "
        "-R {params.genome} "
        "-I {input.tumoral} "
        "-I {input.normal} "
        "-tumor {params.tumor_bam} "
        "-normal {params.normal_bam} "
        "--germline-resource {params.germline_resource} "
        "--af-of-alleles-not-in-resource 0.0000025 "
        "{params.param} "
        "-L {params.intervals} "
        "-O {output.vcf} "
        "-bamout {output.bam} "
        "--native-pair-hmm-threads {threads} "
        "--max-reads-per-alignment-start 0 "
        "--f1r2-tar-gz {output.fir} "
        "--tmp-dir {resources.tmpdir} "
        ">& {log} "


rule learn_orientation_model:
    input:
        rules.mutect_matched.output.fir,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.read-orientation-model.tar.gz",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        exac=config.get("params").get("gatk").get("exac"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/{sample}.pileupsummaries_T.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/mutect/{sample}.pileupsummaries_T.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk LearnReadOrientationModel "
        "--java-options {params.custom} "
        "-I {input} "
        "-O {output} "
        ">& {log} "


rule pileup_summaries_tumoral:
    input:
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.getpileupsummaries.table",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        intervals=config.get("resources").get("bed"),
        exac=config.get("params").get("gatk").get("exac"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/{sample}.pileupsummaries_T.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/mutect/{sample}.pileupsummaries_T.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk GetPileupSummaries "
        "--java-options {params.custom} "
        "-I {input.tumoral} "
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "


rule pileup_summaries_normal:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.normal_pileups.table",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        intervals=config.get("resources").get("bed"),
        exac=config.get("params").get("gatk").get("exac"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/{sample}.pileupsummaries_C.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/mutect/{sample}.pileupsummaries_C.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk GetPileupSummaries "
        "--java-options {params.custom} "
        "-I {input.normal} "
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "


rule calculate_contamination:
    input:
        tab_t=rules.pileup_summaries_tumoral.output,
        tab_c=rules.pileup_summaries_normal.output,
    output:
        table=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.contamination.table",
        ),
        segment=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.tumor.segment",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/{sample}.calculate_contamination.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/mutect/{sample}.calculate_contamination.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk CalculateContamination "
        "--java-options {params.custom} "
        "-I {input.tab_t} "
        "-O {output.table} "
        "-matched {input.tab_c} "
        "--tumor-segmentation {output.segment} "
        ">& {log} "


rule filter_mutect:
    input:
        vcf=rules.mutect_matched.output.vcf,
        cont_tab=rules.calculate_contamination.output.table,
        orient=rules.learn_orientation_model.output,
        segments=rules.calculate_contamination.output.segment,
        stats=rules.mutect_matched.output.stats,
    output:
        vcf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/matched/{sample}.somatic.filtered.vcf.gz",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/{sample}.filter_info.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/mutect/{sample}.filter_info.txt",
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

rule gatk_SelectVariants:
    input:
        lambda wildcards: select_filtered(wildcards),
    output:
        vcf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/mutect/{sample}.somatic.mutect.vcf.gz",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        genome=config.get("resources").get("reference"),
        arguments=config.get("params").get("gatk").get("SelectVariants"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/mutect/{sample}.SelectVariants.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/mutect/{sample}.SelectVariants.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk SelectVariants "
        "--java-options {params.custom} "
        "-R {params.genome} "
        "-V {input} "
        "-O {output.vcf} "
        "{params.arguments} "
        ">& {log} "

rule mutect_hold_on:
    input:
        snvs=rules.gatk_SelectVariants.output.vcf,

    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.mutect.vcf.gz",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "cp {input.snvs} {output.snvs} "

rule mutect_tbi:
    input:
        snvs=rules.mutect_hold_on.output.snvs,
    output:
        snvs = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.mutect.vcf.gz.tbi",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/tabix.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "tabix -p vcf {input.snvs} "
