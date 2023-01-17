rule strelka:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/strelka/{sample}.somatic.strelka.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/strelka/{sample}.somatic.strelka.indels.vcf.gz",
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        normal_name=lambda wildcards,input: get_name(input.normal_name),
        tumor_name=lambda wildcards, input: get_name(input.tumor_name),
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/strelka/{sample}",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/strelka/{sample}.strelka.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/strelka.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input.normal} " 
        "--tumorBam {input.tumoral} "
        "--referenceFasta {params.genome} "
        "--outputCallableRegions --targeted "
        "--callRegions {params.intervals} "
        "--runDir {params.out} ; "
        "python {params.out}/runWorkflow.py -m local -g 10 >& {log} ; "
        "gunzip -c {params.out}/results/variants/somatic.snvs.vcf.gz | sed 's/NORMAL/{params.normal_name}/' | sed 's/TUMOR/{params.tumor_name}/' | bgzip > {output.snvs} ; "
        "gunzip -c {params.out}/results/variants/somatic.indels.vcf.gz | sed 's/NORMAL/{params.normal_name}/' | sed 's/TUMOR/{params.tumor_name}/' | bgzip > {output.indels} "

        #"cp {params.out}/results/variants/somatic.snvs.vcf.gz {output.snvs} && "
        #"cp {params.out}/results/variants/somatic.indels.vcf.gz {output.indels} "
