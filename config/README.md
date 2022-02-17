
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.yaml`. For each sample, specify a `normal_bam` [optional] and one (or more SOON) `tumor_bam`. 
* If `normal_bam` is not provided, the pipeline performs a **tumor-only** analysis.
* You don't have to rename your bam files, the pipeline detect the correct sample name from the BAM header.
