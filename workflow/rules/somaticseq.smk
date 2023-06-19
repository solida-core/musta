rule somaticseq:
    input: lambda wildcards: get_input_files(wildcards),
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
        """
        somaticseq_parallel.py 
        -outdir {params.outdir} 
        -ref {params.genome} 
        #"-dbsnp {params.dbsnp} 
        # "-cosmic "
        --inclusion-region {input.intervals} 
        paired "
        --tumor-bam-file {input.tumoral} 
        --normal-bam-file {input.normal} 
        
        """

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




