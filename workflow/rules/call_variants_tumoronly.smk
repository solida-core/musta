# TUMOR-ONLY CALL WORKFLOW
rule get_tumoronly_sample_names:
    input:
        tumoral= lambda wildcards: get_tumoral_bam(wildcards)
    output:
        tumor=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/tmp/tumoronly/{sample}_tumor.samplename.txt")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"),multiply_by=5)
    log:
        resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"logs/gatk/getsamplename/{sample}.gsn.log")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    shell:
        "gatk "
        "--java-options {params.custom} "
        "GetSampleName "
        "-I {input.tumoral} "
        "-O {output.tumor} "
        ">& {log} "

rule mutect_tumoronly:
    input:
        tumoral= lambda wildcards: get_tumoral_bam(wildcards),
        tumor_name=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/tmp/tumoronly/{sample}_tumor.samplename.txt")
    output:
        vcf=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/tumoronly/{sample}_somatic.vcf.gz"),
        bam=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/tumoronly/{sample}_tumor_normal.bam"),
        fir=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/tumoronly/{sample}_tumor_normal_f1r2.tar.gz"),
        stats=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/tumoronly/{sample}_somatic.vcf.gz.stats")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        genome=config.get("resources").get("reference"),
        intervals=resolve_single_filepath(config.get("paths").get("workdir"),resolve_single_filepath("resources",config.get("resources").get("bed"))),
        param=config.get("params").get("gatk").get("Mutect"),
        germline_resource=config.get("params").get("gatk").get("germline"),
        tumor_bam= lambda wildcards,input: get_name(input.tumor_name)
    log:
        resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"logs/gatk/Mutect2/{sample}.mutect.log")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    shell:
        "gatk "
        "--java-options {params.custom} "
        "Mutect2 "
        "-R {params.genome} "
        "-I {input.tumoral} "
        "--tumor {params.tumor_bam} "
        "--germline-resource {params.germline_resource} "
        "--af-of-alleles-not-in-resource 0.0000025 "
        "{params.param} "
        "-L {params.intervals} "
        "-O {output.vcf} "
        "-bamout {output.bam} "
        "--max-reads-per-alignment-start 0 "
        "--genotype-germline-sites "
        "--f1r2-tar-gz {output.fir} "
        ">& {log} "


rule orientation_model_tumoronly:
    input:
        rules.mutect_tumoronly.output.fir
    output:
        resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/filters/tumoronly/{sample}_read-orientation-model.tar.gz")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        exac=config.get("params").get("gatk").get("exac")
    log:
        resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"logs/gatk/Mutect2/{sample}_orientation.log")
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    shell:
        "gatk LearnReadOrientationModel "
        "--java-options {params.custom} "
        "-I {input} "
        "-O {output} "
        ">& {log} "


rule pileup_summaries_tumoronly:
    input:
        tumoral= lambda wildcards: get_tumoral_bam(wildcards)
    output:
#        resolve_single_filepath(config.get("paths").get("workdir"),"results/filters/tumoronly/{sample}_getpileupsummaries.table")
        resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/filters/tumoronly/{sample}_getpileupsummaries.table")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        intervals=resolve_single_filepath(config.get("paths").get("workdir"),resolve_single_filepath("resources",config.get("resources").get("bed"))),
        exac=config.get("params").get("gatk").get("exac")
    log:
        resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"logs/gatk/Mutect2/{sample}_pileupsummaries_T.log")
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    shell:
        "gatk GetPileupSummaries "
        "--java-options {params.custom} "
        "-I {input.tumoral} "
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "


rule calculate_contamination_tumoronly:
    input:
        tab_t=rules.pileup_summaries_tumoronly.output
    output:
        table=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/filters/tumoronly/{sample}_contamination.table"),
        segment=resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"results/filters/tumoronly/{sample}_tumor.segment")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"),multiply_by=5)
    log:
        resolve_results_filepath(config.get("paths").get("workdir"),config.get("paths").get("project_name"),"logs/gatk/Mutect2/{sample}_calculatecontamination.log")
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    shell:
        "gatk CalculateContamination "
        "--java-options {params.custom} "
        "-I {input.tab_t} "
        "-O {output.table} "
        "--tumor-segmentation {output.segment} "
        ">& {log} "
