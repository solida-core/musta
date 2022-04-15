rule Funcotator:
    input:
        vcf=get_annotation_input()
    output:
        vcf=report(
            resolve_single_filepath(config.get("paths").get("workdir"),"results/annotation/funcotator/{sample}_funcotated.maf"),
            caption="../report/vcf.rst",
            category="Annotation",
        )
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"),multiply_by=5),
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        resources=config.get("params").get("gatk").get("Funcotator")
    log:
        resolve_single_filepath(config.get("paths").get("workdir"),"logs/gatk/Funcotator/{sample}.funcotator.log")
    conda:
       "../envs/gatk.yaml"
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
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
