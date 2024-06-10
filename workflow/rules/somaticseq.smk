rule somaticseq:
    input:
        normal = lambda wildcards: get_normal_bam(wildcards),
        tumoral = lambda wildcards: get_tumoral_bam(wildcards),
        intervals=rules.prepare_bedfile.output.intervals,
        muse=rules.muse_hold_on.output.snvs if config["callers"]["muse"] else [],
        vardict=rules.vardict_hold_on.output.snvs if config["callers"]["vardict"] else [],
        strelka_snvs=rules.strelka_hold_on.output.snvs if config["callers"]["strelka"] else [],
        strelka_indels=rules.strelka_hold_on.output.indels if config["callers"]["strelka"] else [],
        mutect=rules.mutect_hold_on.output.snvs if config["callers"]["mutect"] else [],
        varscan_snvs=rules.varscan_hold_on.output.snvs if config["callers"]["varscan"] else [],
        varscan_indels=rules.varscan_hold_on.output.indels if config["callers"]["varscan"] else [],
        lofreq_snvs=rules.lofreq_hold_on.output.snvs if config["callers"]["lofreq"] else [],
        lofreq_indels=rules.lofreq_hold_on.output.indels if config["callers"]["lofreq"] else [],
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
        mutect_flag=lambda wildcards, input: generate_flag("mutect-vcf", input.mutect),
        vardict_flag=lambda wildcards, input: generate_flag('vardict-vcf', input.vardict),
        muse_flag=lambda wildcards, input: generate_flag('muse-vcf', input.muse),
        varscan_snvs_flag=lambda wildcards, input: generate_flag('varscan-snv', input.varscan_snvs),
        varscan_indels_flag=lambda wildcards, input: generate_flag('varscan-indel', input.varscan_indels),
        lofreq_snvs_flag=lambda wildcards, input: generate_flag('lofreq-snv', input.lofreq_snvs),
        lofreq_indels_flag=lambda wildcards, input: generate_flag('lofreq-indel', input.lofreq_indels),
        strelka_snvs_flag=lambda wildcards, input: generate_flag('strelka-snv', input.strelka_snvs),
        strelka_indels_flag=lambda wildcards, input: generate_flag('strelka-indel', input.strelka_indels),
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
        #"-dbsnp {params.dbsnp} "
        # "-cosmic "
        "--inclusion-region {input.intervals} "
        "paired "
        "--tumor-bam-file {input.tumoral} "
        "--normal-bam-file {input.normal} "
        "{params.mutect_flag} "
        "{params.muse_flag} "
        "{params.vardict_flag} "
        "{params.lofreq_snvs_flag} "
        "{params.lofreq_indels_flag} "
        "{params.strelka_snvs_flag} "
        "{params.strelka_indels_flag} "
        "{params.varscan_snvs_flag} "
        "{params.varscan_indels_flag} "


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
    conda: resolve_single_filepath(config.get("paths").get("workdir"), "workflow/envs/tabix.yaml"),
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
            "detection/results/{sample}.musta.snvs.vcf.gz",
        ),
        indels=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.musta.indels.vcf.gz",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "cp {input.snvs} {output.snvs} ; "
        "cp {input.indels} {output.indels} "

rule somaticseq_tbi:
    input:
        snvs=rules.somaticseq_hold_on.output.snvs,
        indels=rules.somaticseq_hold_on.output.indels,
    output:
        snvs = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.musta.snvs.vcf.gz.tbi",
        ),
        indels = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.musta.indels.vcf.gz.tbi",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/tabix.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "tabix -p vcf {input.snvs} ; "
        "tabix -p vcf {input.indels}  "
