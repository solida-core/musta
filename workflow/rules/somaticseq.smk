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
        mutect="--mutect2-vcf {input.mutect} " if input.mutect else "",
        vardict="--vardict-vcf {input.vardict} " if input.vardict else "",
        muse="--muse-vcf {input.muse} " if input.muse else "",
        varscan_snvs="--varscan-snv {input.varscan_snvs} " if input.varscan_snvs else "",
        varscan_indels="--varscan-indel {input.varscan_indels} " if input.varscan_indels else "",
        lofreq_snvs="--lofreq-snv {input.lofreq_snvs} " if input.lofreq_snvs else "",
        lofreq_indels="--lofreq-indel {input.lofreq_indels} " if input.lofreq_indels else "",
        strelka_snvs="--strelka-snv {input.strelka_snvs} " if input.strelka_snvs else "",
        strelka_indels="--strelka-indel {input.strelka_indels} " if input.strelka_indels else "",
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

rule somaticseq_tbi:
    input:
        snvs=rules.somaticseq_hold_on.output.snvs,
        indels=rules.somaticseq_hold_on.output.indels,
    output:
        snvs = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.consensus.snvs.vcf.gz.tbi",
        ),
        indels = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.consensus.indels.vcf.gz.tbi",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/tabix.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "tabix -p vcf {input.snvs} ; "
        "tabix -p vcf {input.indels}  "
