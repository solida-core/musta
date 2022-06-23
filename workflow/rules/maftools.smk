
rule maftools:
    input:
        mafs=get_maf_file_input()
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/annotation/maftools/signature_contributions.png")
    params:
        project_id="prova",
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/annotation/maftools")
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/scripts/maftools.R")


rule maftools_general:
    input:
        mafs=get_maf_file_input()
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/general/plots/top10_VAF.png")
    params:
        project_id="prova",
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/general")
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
        mafs=get_maf_file_input()
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/signatures/plots/cosmic_signatures.png")
    params:
        project_id="prova",
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/signatures")
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
        mafs=get_maf_file_input()
    output:
        resolve_single_filepath(config.get("paths").get("results_dir"),"results/driver/plots/somatic_interactions.png")
    params:
        outdir=resolve_single_filepath(config.get("paths").get("results_dir"),"results/driver")
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/scripts/maftools_driver.R")
