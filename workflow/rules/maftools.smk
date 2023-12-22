rule maftools_base:
    input:
        mafs=get_maf_file_input(),
    output:
        summary=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/base/plots/summary.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/summary.rst"
            ),
            category="Maftools",
        ),
        oncoplot=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/base/plots/large_oncoplot.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/oncoplot.rst"
            ),
            category="Maftools",
        ),
        tgca=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/base/plots/TGCA_compare.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/tgca.rst"
            ),
            category="Maftools",
        ),
        vaf=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/base/plots/top10_VAF.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/vaf.rst"
            ),
            category="Maftools",
        ),
    params:
        project_id=config.get("params").get("maftools").get("project_name"),
        outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/base.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/base.txt",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/maftools_basic.R"
        )


rule maftools_signatures:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary,
    output:
        signatures=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/signatures/plots/cosmic_signatures.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/signatures.rst"
            ),
            category="Signatures",
        ),
    params:
        project_id=config.get("params").get("maftools").get("project_name"),
        outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/signatures.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/signatures.txt",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/maftools_signatures.R"
        )


rule maftools_driver:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary,
    output:
        interactions=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/driver/plots/somatic_interactions.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/interactions.rst"
            ),
            category="Somatic Interactions",
        ),
    params:
        outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/gene_drivers.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/gene_drivers.txt",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/maftools_driver.R"
        )


rule maftools_pathways:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary,
    output:
        resolve_single_filepath(
            config.get("paths").get("results_dir"),
            "interpretation/pathways/plots/oncogenic_pathways.png",
        ),
    params:
        outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/pathways.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/pathways.txt",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/maftools_pathway.R"
        )


rule maftools_heterogeneity:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary,
    output:
        resolve_single_filepath(
            config.get("paths").get("results_dir"),
            "interpretation/heterogeneity/tables/successful.tsv",
        ),
    params:
        outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/heterogeneity.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/heterogeneity.txt",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"),
            "workflow/scripts/maftools_heterogeneity.R",
        )
