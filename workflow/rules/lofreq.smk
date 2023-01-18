rule prepare_bedfile:
    input:
        intervals=config.get("resources").get("bed"),
    output:
        intervals=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/bedfile/bed.vcf",
        ),
    shell:
        "gunzip -c {input.intervals} > {output.intervals}"

rule lofreq:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/lofreq/{sample}.somatic.lofreq.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/lofreq/{sample}.somatic.lofreq.indels.vcf.gz",
        ),
    params:
        genome=config.get("resources").get("reference"),
        intervals=rules.prepare_bedfile.output.intervals,
        dbsnp=config.get("resources").get("dbsnp"),
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "variant_calling/lofreq/{sample}.",
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/variant_calling/lofreq/{sample}.lofreq.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/lofreq.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "lofreq somatic "
        "-f {params.genome} "
        "-t {input.tumoral} "
        "-n {input.normal} "
        "--threads {threads} "
        "-l {params.intervals} "
        "-d {params.dbsnp} "
        "-o {params.out} "
        "--call-indels "
        ">& {log} ; "
        "cp {params.out}somatic_final_minus-dbsnp.snvs.vcf.gz {output.snvs} && "
        "cp {params.out}somatic_final_minus-dbsnp.indels.vcf.gz {output.indels} "

