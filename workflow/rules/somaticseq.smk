rule somaticseq:
    input:
        normal = lambda wildcards: get_normal_bam(wildcards),
        tumoral = lambda wildcards: get_tumoral_bam(wildcards),
        muse=muse_flag,
        vardict=vardict_flag,
        strelka_snvs=strelka_snvs_flag,
        strelka_indels=strelka_indels_flag,
        mutect=mutect_flag,
        varscan_snvs=varscan_snvs_flag,
        varscan_indels=varscan_indels_flag,
        lofreq_snvs=lofreq_snvs_flag,
        lofreq_indels=lofreq_indels_flag,
        intervals=rules.prepare_bedfile.output.intervals,
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
        "{input.mutect} "
        "{input.varscan_snvs} "
        "{input.varscan_indels} "
        "{input.vardict} "
        "{input.muse} "
        "{input.lofreq_snvs} "
        "{input.lofreq_indels} "
        "{input.strelka_snvs} "
        "{input.strelka_indels} "

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

#############################################
def mutect_flag(wildcards):
    if "mutect_hold_on" in rules:
        return "--mutect2-vcf {}".format(rules.mutect_hold_on.output.snvs)
    else:
        return ""

def varscan_snvs_flag(wildcards):
    if "varscan_hold_on" in rules:
        return "--varscan-snv {}".format(rules.varscan_hold_on.output.snvs)
    else:
        return ""

# Definisci le altre funzioni per gli input condizionali

# Nella definizione della regola, aggiungi i flag agli input corrispondenti

rule somaticseq:
    input:
        normal = lambda wildcards: get_normal_bam(wildcards),
        tumoral = lambda wildcards: get_tumoral_bam(wildcards),
        intervals=rules.prepare_bedfile.output.intervals,
        mutect_flag = mutect_flag,
        varscan_snvs_flag = varscan_snvs_flag,
        # Aggiungi gli altri flag
    output:
        # Specifica gli output della regola
    params:
        outdir = "path/to/outdir",
        genome = "path/to/genome",
        dbsnp = "path/to/dbsnp"
    shell:
        """
        somaticseq_parallel.py
        -outdir {params
