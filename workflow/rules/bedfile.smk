rule prepare_bedfile:
    input:
        intervals=config.get("resources").get("bed"),
    output:
        intervals=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "bedfile/bedfile.vcf",
        ),
    shell:
        "gunzip -c {input.intervals} > {output.intervals}"
