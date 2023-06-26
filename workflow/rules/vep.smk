rule vep:
    input:
        vcf=lambda wildcards: get_annotation_input(wildcards),
        normal_name= rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        vcf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/vep/{sample}.annotated.vep.vcf",
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
            "classification/vep/{sample}.vep.log",
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
        "vep "
        "--input_file {input.vcf} "
        "--output_file {output.vcf} "
        "--fasta {params.genome} "
        "--dir {params.resources} "
        "--assembly {params.genome_version} "
        "--cache --cache_version {params.cache_version} "
        "--sift b --ccds --uniprot --hgvs --symbol "
        "--numbers --domains --gene_phenotype --canonical "
        "--protein --biotype --tsl --pubmed --variant_class "
        "--shift_hgvs 1 --check_existing --total_length "
        "--allele_number --no_escape --xref_refseq --failed 1 "
        "--vcf "
        "--minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length " 
        "--polyphen b --af --af_1kg --af_esp --regulatory "
        ">& {log} "

rule vep2maf:
    input:
        vcf=rules.vep.output.vcf,
        normal_name= rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        maf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/vep/{sample}.annotated.vep.maf",
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
            "classification/vep/{sample}.vep2maf.log",
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


rule vep_hold_on:
    input:
        vcf=rules.vep.output.vcf,
        maf=rules.vep2maf.output.maf,
    output:
        vcf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/results/{sample}.annotated.vep.vcf.gz",
        ),
        tbi=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/results/{sample}.annotated.vep.vcf.gz.tbi",
        ),
        maf=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "classification/results/{sample}.annotated.vep.maf",
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


