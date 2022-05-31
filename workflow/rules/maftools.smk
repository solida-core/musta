
rule maftools:
    input:
        expand(resolve_results_filepath(config.get("paths").get("results_dir"),"results/annotation/funcotator/{sample}_funcotated.maf"),sample=list(samples_master.keys()))
    output:
        resolve_single_filepath(config.get("paths").get("workdir"),"results/annotation/maftools/signature_contributions.png")
    params:
        project_id="prova",
        outdir=resolve_single_filepath(config.get("paths").get("workdir"),"results/annotation/maftools")
    conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/maftools.yaml")
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    script:
        "scripts/maftools.R"
