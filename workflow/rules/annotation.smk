rule Funcotator:
    input:
        vcf="results/{sample}_somatic_filtered_selected.vcf.gz"
    output:
        vcf="results/annotation/funcotator/{sample}_funcotated.maf"
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"),multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path("ref"),config.get("ref").get("fasta")),
        intervals=config.get("processing").get("interval_list"),
        resources=config.get("params").get("gatk").get("Funcotator")
    log:
        "logs/gatk/Funcotator/{sample}.funcotator.log"
    conda:
       "../envs/gatk.yaml"
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk Funcotator "
        "--java-options {params.custom} "
        "-R {params.genome} "
        "-L {params.intervals} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--data-sources-path {params.resources} "
        "--output-file-format MAF "
        "--ref-version hg19 "
        ">& {log} "


rule kggseq:
    input:
        vcf="results/{sample}_somatic_filtered_selected.vcf.gz"
    output:
        vcf="results/annotation/kggseq/{sample}.selected.flt.vcf",
        log="results/annotation/kggseq/{sample}selected.log",
        txt="results/annotation/kggseq/{sample}selected.flt.txt"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        software=config.get("params").get("kggseq").get("software"),
        params=config.get("params").get("kggseq").get("parameters"),
        out_basename= lambda wildcards,output: output.vcf[:-8]
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "{params.software} "
        "-nt {threads} "
        "{params.custom} "
        "{params.params} "
        "--no-gty-vcf-file {input.vcf} "
        "--out {params.out_basename} "
