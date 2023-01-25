rule somaticseq:
    input:
        normal = lambda wildcards: get_normal_bam(wildcards),
        tumoral = lambda wildcards: get_tumoral_bam(wildcards),
        muse=rules.muse_hold_on.output.snvs,
        vardict=rules.vardict_hold_on.output.snvs,
        strelka_snvs=rules.strelka_hold_on.output.snvs,
        strelka_indels=rules.strelka_hold_on.output.indels,
        mutect=rules.mutect_hold_on.output.snvs,
        varscan_snvs=rules.varscan_hold_on.output.snvs,
        varscan_indels=rules.varscan_hold_on.output.indels,
        lofreq_snvs=rules.lofreq_hold_on.output.snvs,
        lofreq_indels=rules.lofreq_hold_on.output.indels
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/somaticseq/{sample}/Consensus.sSNV.vcf",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/somaticseq/{sample}/Consensus.sINDEL.vcf",
        ),
    params:
        outdir=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/somaticseq/{sample}",
        ),
        genome=config.get("resources").get("reference"),
        dbsnp=config.get("resources").get("dbsnp"),
        intervals=config.get("resources").get("bed"),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "detection/somaticseq/{sample}.somaticseq.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/somatiqseq/{sample}.somaticseq.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/somaticseq.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "somaticseq_parallel.py "
        "-outdir {params.outdir} "
        "-ref {params.genome} "
        "-dbsnp {params.dbsnp} "
        # "-cosmic "
        "--inclusion-region {params.intervals} "
        "paired "
        "--tumor-bam-file {input.tumoral} "
        "--normal-bam-file {input.normal} "
        "--mutect2-vcf {input.mutect} "
        "--varscan-snv {input.varscan_snvs} "
        "--varscan-indel {input.varscan_indels} "
        "--vardict-vcf {input.vardict} "
        "--muse-vcf {input.muse} "
        "--lofreq-snv {input.lofreq_snvs} "
        "--lofreq-indel {input.lofreq_indels} "
        "--strelka-snv {input.strelka_snvs} "
        "--strelka-indel {input.strelka_indels} "

rule somaticseq_out:
    input:
        snvs=rules.somaticseq.output.snvs,
        indels=rules.somaticseq.output.indels,
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/somaticseq/{sample}/{sample}.somaticseq.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/somaticseq/{sample}/{sample}.somaticseq.indels.vcf.gz",
        ),
    conda: resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "bgzip -c {input.snvs} > {output.snvs} ; "
        "bgzip -c {input.indels} > {output.indels} "

rule somaticseq_hold_on:
    input:
        snvs=rules.somaticseq_out.output.snvs,
        indels=rules.somaticseq_out.output.indels,
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.consensus.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.consensus.indels.vcf.gz",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "cp {input.snvs} {output.snvs} ; "
        "cp {input.indels} {output.indels} "
