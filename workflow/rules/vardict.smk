rule vardict:
    input:
        normal=lambda wildcards: get_normal_bam(wildcards),
        tumoral=lambda wildcards: get_tumoral_bam(wildcards),
        normal_name=rules.get_sample_names.output.normal,
        tumor_name=rules.get_sample_names.output.tumor,
        intervals=rules.prepare_bedfile.output.intervals,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/vardict/{sample}.somatic.vardict.snvs.vcf.gz",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/vardict.yaml"
        ),
    params:
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/vardict/{sample}.somatic.vardict.vcf",
        ),
        genome=config.get("resources").get("reference"),
        normal_name= lambda wildcards,input: get_name(input.normal_name),
        tumor_name=lambda wildcards, input: get_name(input.tumor_name),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "detection/vardict/{sample}.vardict.txt",
        ),
    shell:
        r"""
        set -euo pipefail

        # Locate VarDict auxiliary scripts
        VARDICT_DIR="$(dirname "$(readlink -f "$(command -v vardict)")")"

        test -f "$VARDICT_DIR/testsomatic.R"
        test -f "$VARDICT_DIR/var2vcf_paired.pl"

        # Somatic variant calling (tumor-normal)
        vardict \
            -G {params.genome:q} \
            -f 0.01 \
            -N {params.tumor_name:q} \
            -b "{input.tumoral}|{input.normal}" \
            -c 1 -S 2 -E 3 \
            {input.intervals:q} \
        | Rscript "$VARDICT_DIR/testsomatic.R" \
        | perl "$VARDICT_DIR/var2vcf_paired.pl" \
            -N "{params.tumor_name}|{params.normal_name}" \
            -f 0.01 \
        > {params.out:q}

        bgzip -c {params.out:q} > {output:q}
        """

rule vardict_hold_on:
    input:
        snvs=rules.vardict.output,
    output:
        snvs=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.vardict.snvs.vcf.gz"
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "cp {input.snvs} {output.snvs} "

rule vardict_tbi:
    input:
        snvs=rules.vardict_hold_on.output.snvs,
    output:
        snvs = resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "detection/results/{sample}.somatic.vardict.snvs.vcf.gz.tbi",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/tabix.yaml"
        ),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99),
    shell:
        "tabix -p vcf {input.snvs} "
