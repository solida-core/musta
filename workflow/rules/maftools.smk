
rule maftools_base:
    input:
        mafs=get_maf_file_input()
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/general/plots/top10_VAF.png")
    params:
        project_id="prova",
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/general")
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
        requisites=rules.maftools_base.output
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/signatures/plots/cosmic_signatures.png")
    params:
        project_id="prova",
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
        requisites=rules.maftools_base.output
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/analysis/driver/plots/somatic_interactions.png")
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
        requisites=rules.maftools_base.output
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
        requisites=rules.maftools_base.output
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
