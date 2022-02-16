
### mettere un if con attivazione in caso di matched
rule get_sample_names:
    input:
        lambda wildcards: get_bams(wildcards)
    output:
        normal="results/tmp/{sample}_normal.samplename.txt",
        tumor="results/tmp/{sample}_tumor.samplename.txt",
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"),multiply_by=5)
    log:
        normal="logs/gatk/getsamplename/{sample}.gsn.log",
        tumor="logs/gatk/getsamplename/{sample}.gsn_tumor.log"
    conda:
        "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    shell:
        "gatk "
        "--java-options {params.custom} "
        "GetSampleName "
        "-I {input[1]} " # get tumor sample bam
        "-O {output.tumor} "
        ">& {log.tumor} && "
        "gatk "
        "--java-options {params.custom} "
        "GetSampleName "
        "-I {input[0]} "  # get normal sample bam
        "-O {output.normal} "
        ">& {log.normal} "

rule mutect_matched:
    input:
        lambda wildcards: get_bams(wildcards),
        normal_name="results/tmp/{sample}_normal.samplename.txt",
        tumor_name="results/tmp/{sample}_tumor.samplename.txt"
    output:
        vcf="results/matched/{sample}_somatic.vcf.gz",
        bam="results/matched/{sample}_tumor_normal.bam",
        fir="results/matched/{sample}_tumor_normal_f1r2.tar.gz"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path("ref"), config.get("ref").get("fasta")),
        intervals=config.get("processing").get("interval_list"),
        param=config.get("params").get("gatk").get("Mutect"),
        germline_resource=config.get("germline"),
        normal_bam = ""#funzione_apri_e_leggi("results/tmp/{sample}_normal.samplename.txt")
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
        "-I {input[1]} "#tumoral bam
        "-I {input[0]} "#control bam
        "-normal {params.normal_bam} "#control name
#        "-pon {params.pon} "
        "--germline-resource {params.germline_resource} "
        "--af-of-alleles-not-in-resource 0.0000025 "
        "{params.param} "
        "-L {params.intervals} "
        "-O {output.vcf} "
        "-bamout {output.bam} "
        "--native-pair-hmm-threads {threads} "
        "--max-reads-per-alignment-start 0 "
        "--f1r2-tar-gz {output.fir} "
        ">& {log} "

rule learn_orientation_model:
    input:
        rules.mutect_matched.output.fir
    output:
        "results/filters/matched/{sample}_read-orientation-model.tar.gz"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        exac=config.get("exac")
    log:
        "logs/gatk/Mutect2/{sample}_pileupsummaries_T.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk LearnReadOrientationModel "
        "--java-options {params.custom} "
        "-I {input} "
        "-O {output} "
        ">& {log} "

rule pileup_summaries_tumoral:
    input:
        lambda wildcards: get_bams(wildcards)
    output:
        "results/filters/matched/{sample}_getpileupsummaries.table"
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
        "-I {input[1]} " # corresponding to tbam
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "

rule pileup_summaries_normal:
    input:
        lambda wildcards: get_bams(wildcards)
    output:
        "results/filters/matched/{sample}_normal_pileups.table"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"), multiply_by=5),
        intervals=config.get("processing").get("interval_list"),
        exac=config.get("exac")
    log:
        "logs/gatk/Mutect2/{sample}_pileupsummaries_C.log"
    conda:
       "../envs/gatk.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk GetPileupSummaries "
        "--java-options {params.custom} "
        "-I {input[0]} "
        "-V {params.exac} "
        "-L {params.intervals} "
        "-O {output} "
        ">& {log} "

rule calculate_contamination:
    input:
        tab_t=rules.pileup_summaries_tumoral.output,
        tab_c=rules.pileup_summaries_normal.output
    output:
        table="results/filters/matched/{sample}_contamination.table",
        segment="results/filters/matched/{sample}_tumor.segment"
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
        "-matched {input.tab_c} "
        "--tumor-segmentation {output.segment} "
        ">& {log} "
