rule Funcotator:
    input:
        vcf=get_annotation_input(),
        normal_name= rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        vcf=report(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "classification/results/{sample}.annotated.vcf",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/vcf.rst"
            ),
            category="Annotation",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        resources=config.get("params").get("gatk").get("Funcotator").get("resources"),
        genome_version=config.get("params")
        .get("gatk")
        .get("Funcotator")
        .get("reference_version"),
        normal_name = lambda wildcards, input: get_name(input.normal_name),
        tumor_name = lambda wildcards, input: get_name(input.tumor_name),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "classification/{sample}.funcotator.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/gatk.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "gatk Funcotator "
        "--java-options {params.custom} "
        "-R {params.genome} "
        #"-L {params.intervals} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--data-sources-path {params.resources} "
        "--output-file-format VCF "
        "--annotation-default normal_barcode:{params.normal_name} "
        "--annotation-default tumor_barcode:{params.tumor_name} "
        "--ref-version {params.genome_version} "
        "--cache-version 106 "
        "--tmp-dir {resources.tmpdir} "
        ">& {log} "

rule funcotator_vcf2maf:
    input:
        vcf=rules.Funcotator.output.vcf,
    output:
        maf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/results/{sample}.annotated.maf",
        ),
    params:
        genome=config.get("resources").get("reference"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "classification/{sample}.vcf2maf.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/vcf2maf.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "VEPPATH=$(dirname $(which vep)) ; "
        "vcf2maf.pl "
        "--input-vcf {input.vcf} "
        "--output-maf {output.maf} "
        "--ref-fasta {params.genome} "
        "--tmp-dir {resources.tmpdir} "
        "--vep-path $VEPPATH "
        ">& {log} "
