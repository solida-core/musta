rule funcotator_vcf:
    input:
        vcf=lambda wildcards: get_annotation_input(wildcards),
        normal_name= rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        vcf=report(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "classification/funcotator/{sample}.annotated.funcotator.vcf",
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
        genome_version=config.get("params").get("gatk").get("Funcotator").get("reference_version"),
        normal_name = lambda wildcards, input: get_name(input.normal_name),
        tumor_name = lambda wildcards, input: get_name(input.tumor_name),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "classification/funcotator/{sample}.funcotator.vcf.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "classification/funcotator/{sample}.funcotator.vcf.txt",
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
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--data-sources-path {params.resources} "
        "--output-file-format VCF "
        "--annotation-default normal_barcode:{params.normal_name} "
        "--annotation-default tumor_barcode:{params.tumor_name} "
        "--ref-version {params.genome_version} "
        "--tmp-dir {resources.tmpdir} "
        ">& {log} "


rule funcotator_maf:
    input:
        vcf=rules.funcotator_vcf.output.vcf,
        normal_name= rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        maf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/funcotator/{sample}.annotated.funcotator.maf",
        ),
    params:
        genome=config.get("resources").get("reference"),
        genome_version=get_vep_genome_version(config.get("params").get("vep").get("reference_version")),
        resources=config.get("params").get('vep').get("resources"),
        cache_version=config.get("params").get('vep').get("cache_version"),
        normal_name= lambda wildcards,input: get_name(input.normal_name),
        tumor_name=lambda wildcards, input: get_name(input.tumor_name),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "classification/funcotator/{sample}.funcotator.vcf2maf.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "classification/funcotator/{sample}.funcotator.vcf2maf.txt",
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
        "--normal-id {params.normal_name} "
        "--tumor-id {params.tumor_name} "
        "--ncbi-build {params.genome_version} "
        "--tmp-dir {resources.tmpdir} "
        "--cache-version {params.cache_version} "
        "--vep-path $VEPPATH "
        "--vep-data {params.resources} "
        "--verbose "
        ">& {log} "


rule funcotator_hold_on:
    input:
        vcf=rules.funcotator_vcf.output.vcf,
        maf=rules.funcotator_maf.output.maf,
    output:
        vcf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/results/{sample}.annotated.funcotator.vcf.gz",
        ),
        tbi=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/results/{sample}.annotated.funcotator.vcf.gz.tbi",
        ),
        maf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/results/{sample}.annotated.funcotator.maf",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/tabix.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "cp {input.maf} {output.maf} ; "
        "bgzip -c {input.vcf} > {output.vcf} ; "
        "tabix -p vcf {output.vcf} ; "
