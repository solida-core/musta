## filter tumor-normal
rule filter_mutect:
    input:
        vcf=rules.mutect_matched.output.vcf,
        cont_tab=rules.calculate_contamination.output.table,
        orient=rules.learn_orientation_model.output,
        segments=rules.calculate_contamination.output.segment
    output:
        vcf="results/matched/{sample}_somatic_filtered.vcf.gz"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path("ref"), (config.get("ref").get("fasta"))),
        intervals=config.get("processing").get("interval_list"),
    log:
        "logs/gatk/Mutect2/{sample}.filter_info.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk FilterMutectCalls "
        "--java-options {params.custom} "
        "-V {input.vcf} "
        "--tumor-segmentation {input.segments} "
        "--contamination-table {input.cont_tab} "
        "--ob-priors {input.orient} "
        "-O {output} "
        "-R {params.genome} "
        "-L {params.intervals} "
        ">& {log} "

# filter tumor-only
rule filter_mutect_tumoronly:
    input:
        vcf = rules.mutect_tumoronly.output.vcf,
        cont_tab = rules.calculate_contamination_tumoronly.output.table,
        orient = rules.orientation_model_tumoronly.output,
        segments = rules.calculate_contamination_tumoronly.output.segment
    output:
        vcf="results/tumoronly/{sample}_somatic_filtered.vcf.gz"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path("ref"), (config.get("ref").get("fasta"))),
        intervals=config.get("processing").get("interval_list"),
    log:
        "logs/gatk/Mutect2/{sample}.filter_info.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk FilterMutectCalls "
        "--java-options {params.custom} "
        "-V {input.vcf} "
        "--tumor-segmentation {input.segments} "
        "--contamination-table {input.cont_tab} "
        "--ob-priors {input.orient} "
        "-O {output} "
        "-R {params.genome} "
        "-L {params.intervals} "
        ">& {log} "

def select_filtered(wildcards,samples):
    if samples.loc[wildcards.sample,["normal_bam"]].isna().all():
        # print("tumoronly")
        # TUMOR-ONLY
        return rules.filter_mutect_tumoronly.output.vcf
    # PE
    else:
        return rules.filter_mutect.output.vcf


rule gatk_SelectVariants:
    input:
        # vcf="data/results/{sample}_somatic_filtered.vcf.gz"
        lambda wildcards: select_filtered(wildcards, samples)
    output:
        vcf="results/{sample}_somatic_filtered_selected.vcf.gz"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"),multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path("ref"), (config.get("ref").get("fasta"))),
        arguments=config.get("params").get("gatk").get("SelectVariants")
    log:
        "logs/gatk/SelectVariants/{sample}.SelectVariants.log"
    conda:
       "../envs/gatk.yaml"

    benchmark:
        "benchmarks/gatk/SelectVariants/{sample}.SelectVariants.txt"

    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk SelectVariants "
        "--java-options {params.custom} "
        "-R {params.genome} "
        "-V {input} "
        "-O {output.vcf} "
        "{params.arguments} "
        ">& {log} "
