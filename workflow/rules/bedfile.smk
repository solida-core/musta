rule prepare_bedfile:
    input:
        intervals=config.get("resources").get("bed"),
    output:
        intervals=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "bedfile/bedfile.vcf",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99),
    shell:
        "gunzip -c {input.intervals} > {output.intervals}"
