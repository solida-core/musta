rule strelka:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/strelka/{sample}.somatic.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/strelka/{sample}.somatic.indels.vcf.gz",
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/rollcall/strelka/{sample}",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/rollcall/{sample}.strelka.log",
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
        "python {params.out}/runWorkflow.py -m local -g 10 >& ${log} ; "
        "cp {params.out}/results/variants/somatic.snvs.vcf.gz {output.snvs} && "
        "cp {params.out}/results/variants/somatic.indels.vcf.gz {output.indels} "
