# TUMOR-ONLY CALL WORKFLOW
rule get_tumoronly_sample_names:
    input:
        lambda wildcards: get_bams(wildcards,samples)
    output:
        tumor="results/tmp/tumoronly/{sample}_tumor.samplename.txt",
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"),multiply_by=5)
    log:
        "logs/gatk/getsamplename/{sample}.gsn.log"
    conda:
        "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    shell:
        "gatk "
        "--java-options {params.custom} "
        "GetSampleName "
        "-I {input[0]} " # get tumor sample bam
        "-O {output.tumor} "
        ">& {log} "

rule mutect_tumoronly:
    input:
        lambda wildcards: get_bams(wildcards,samples),
        tumor_name="results/tmp/tumoronly/{sample}_tumor.samplename.txt"
    output:
        vcf="results/tumoronly/{sample}_somatic.vcf.gz",
        bam="results/tumoronly/{sample}_tumor_normal.bam",
        fir="results/tumoronly/{sample}_tumor_normal_f1r2.tar.gz"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path("ref"), config.get("ref").get("fasta")),
        intervals=config.get("processing").get("interval_list"),
        param=config.get("params").get("gatk").get("Mutect"),
        germline_resource=config.get("germline")
    log:
        "logs/gatk/Mutect2/{sample}.mutect.log"
    conda:
        "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk "
        "--java-options {params.custom} "
        "Mutect2 "
        "-R {params.genome} "
        "-I {input[0]} "#tumoral bam
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
        "results/filters/tumoronly/{sample}_read-orientation-model.tar.gz"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        exac=config.get("exac")
    log:
        "logs/gatk/Mutect2/{sample}_orientation.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk LearnReadOrientationModel "
        "--java-options {params.custom} "
        "-I {input} "
        "-O {output} "
        ">& {log} "


rule pileup_summaries_tumoronly:
    input:
        lambda wildcards: get_bams(wildcards,samples)
    output:
        "results/filters/tumoronly/{sample}_getpileupsummaries.table"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        intervals=config.get("processing").get("interval_list"),
        exac=config.get("exac")
    log:
        "logs/gatk/Mutect2/{sample}_pileupsummaries_T.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk GetPileupSummaries "
        "--java-options {params.custom} "
        "-I {input[0]} " # corresponding to tbam
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "


rule calculate_contamination_tumoronly:
    input:
        tab_t=rules.pileup_summaries_tumoronly.output
    output:
        table="results/filters/tumoronly/{sample}_contamination.table",
        segment="results/filters/tumoronly/{sample}_tumor.segment"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"),multiply_by=5)
    log:
        "logs/gatk/Mutect2/{sample}_calculatecontamination.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk CalculateContamination "
        "--java-options {params.custom} "
        "-I {input.tab_t} "
        "-O {output.table} "
        "--tumor-segmentation {output.segment} "
        ">& {log} "
