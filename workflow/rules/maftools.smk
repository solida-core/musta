rule maftools_base:
    input:
        mafs=get_maf_file_input(),
    output:
        summary=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/results/variant_visualization/plots/summary.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/summary.rst"
            ),
            category="Maftools",
        ),
        oncoplot=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/results/variant_visualization/plots/large_oncoplot.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/oncoplot.rst"
            ),
            category="Maftools",
        ),
        tgca=report(
            resolve_single_filepath(
                config.get("paths").get("results_dir"),
                "interpretation/results/variant_visualization/plots/TGCA_compare.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/tgca.rst"
            ),
            category="Maftools",
        ),

    params:
        project_id=config.get("params").get("maftools").get("project_name"),
        all_variants=config.get("variants").get("all"),
        # outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
        outdir=resolve_single_filepath(
            config.get("paths").get("results_dir"),
            "interpretation/results/variant_visualization",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/variant_visualization.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/variant_visualization.txt",
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
                "interpretation/results/mutational_signature_analysis/plots/cosmic_signatures.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/signatures.rst"
            ),
            category="Signatures",
        ),
    params:
        project_id=config.get("params").get("maftools").get("project_name"),
        all_variants=config.get("variants").get("all"),
        # outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
        outdir=resolve_single_filepath(
            config.get("paths").get("results_dir"),
            "interpretation/results/mutational_signature_analysis",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/mutational_signature_analysis.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/mutational_signature_analysis.txt",
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
                "interpretation/results/driver_gene_identification/plots/oncodrive.png",
        ),
        caption=resolve_single_filepath(
        config.get("paths").get("workdir"), "workflow/report/oncodrive.rst"
            ),
            category="Somatic Interactions",
        ),
    params:
        all_variants=config.get("variants").get("all"),
        # outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
        outdir=resolve_single_filepath(
            config.get("paths").get("results_dir"),
            "interpretation/results/driver_gene_identification",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/driver_gene_identification.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/driver_gene_identification.txt",
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
            "interpretation/results/pathway_analysis/plots/oncogenic_pathways.png",
        ),
    params:
        all_variants=config.get("variants").get("all"),
        # outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
        outdir=resolve_single_filepath(
            config.get("paths").get("results_dir"),
            "interpretation/results/pathway_analysis",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/pathway_analysis.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/pathway_analysis.txt",
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
            "interpretation/results/tumor_heterogeneity/tables/successful.tsv",
        ),
    params:
        all_variants=config.get("variants").get("all"),
        # outdir=lambda w, output: os.path.split(os.path.split(output[0])[0])[0],
        outdir=resolve_single_filepath(
            config.get("paths").get("results_dir"),
            "interpretation/results/tumor_heterogeneity",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/maftools.yaml"
        )
    log:
        resolve_results_filepath(
            config.get("paths").get("log_dir"),
            "interpretation/tumor_heterogeneity.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("bench_dir"),
            "interpretation/tumor_heterogeneity.txt",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"),
            "workflow/scripts/maftools_heterogeneity.R",
        )
