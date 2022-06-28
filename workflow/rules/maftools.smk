
rule maftools_base:
    input:
        mafs=get_maf_file_input()
    output:
        summary=report(
            resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/base/plots/summary.png"),
            caption=resolve_single_filepath(config.get("paths").get("workdir"),"workflow/report/summary.rst"),
            category="Maftools"
        ),
        oncoplot=report(
            resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/base/plots/large_oncoplot.png"),
            caption=resolve_single_filepath(config.get("paths").get("workdir"),"workflow/report/oncoplot.rst"),
            category="Maftools"
        ),
        tgca=report(
            resolve_single_filepath(config.get("paths").get("results_dir"), "results/analysis/base/plots/TGCA_compare.png"),
            caption=resolve_single_filepath(config.get("paths").get("workdir"),"workflow/report/tgca.rst"),
            category="Maftools"
        ),
        vaf=report(
            resolve_single_filepath(config.get("paths").get("results_dir"), "results/analysis/base/plots/top10_VAF.png"),
            caption=resolve_single_filepath(config.get("paths").get("workdir"),"workflow/report/vaf.rst"),
            category="Maftools"
        )
    params:
        project_id=config.get("params").get("maftools").get("project_name"),
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/base")
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/scripts/maftools_basic.R")


rule maftools_signatures:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary
    output:
        signatures=report(
            resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/signatures/plots/cosmic_signatures.png"),
            caption=resolve_single_filepath(config.get("paths").get("workdir"),"workflow/report/signatures.rst"),
            category="Signatures",
        )
    params:
        project_id=config.get("params").get("maftools").get("project_name"),
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/signatures")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/scripts/maftools_signatures.R")


rule maftools_driver:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary
    output:
        interactions=report(
            resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/driver/plots/somatic_interactions.png"),
            caption=resolve_single_filepath(config.get("paths").get("workdir"),"workflow/report/interactions.rst"),
            category="Somatic Interactions",
        )
    params:
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/driver")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/scripts/maftools_driver.R")


rule maftools_pathways:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/pathways/plots/oncogenic_pathways.png")
    params:
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/pathways")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/scripts/maftools_pathway.R")


rule maftools_heterogeneity:
    input:
        mafs=get_maf_file_input(),
        requisites=rules.maftools_base.output.summary
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/heterogeneity/tables/successful.tsv")
    params:
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/heterogeneity")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/scripts/maftools_heterogeneity.R")
