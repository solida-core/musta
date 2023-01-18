rule strelka:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/strelka/{sample}/results/{sample}.somatic.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/strelka/{sample}/results/{sample}.somatic.indels.vcf.gz",
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
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
        "python {params.out}/runWorkflow.py -m local -g 10 >& {log}  "

rule strelka_out:
    input:
        snvs=rules.strelka.output.snvs,
        indels=rules.strelka.output.indels,
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
        normal_name=lambda wildcards,input: get_name(input.normal_name),
        tumor_name=lambda wildcards, input: get_name(input.tumor_name),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/samtools.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "gunzip -c {input.snvs} | sed 's/NORMAL/{params.normal_name}/' | sed 's/TUMOR/{params.tumor_name}/' | bgzip -c > {output.snvs} ; "
        "gunzip -c {input.indels} | sed 's/NORMAL/{params.normal_name}/' | sed 's/TUMOR/{params.tumor_name}/' | bgzip -c > {output.indels} "
