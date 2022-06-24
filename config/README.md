
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.yaml`. For each sample, specify a `normal_bam` [optional] and one (or more SOON) `tumor_bam`. 
* If `normal_bam` is not provided, the pipeline performs a **tumor-only** analysis.
* You don't have to rename your bam files, the pipeline detect the correct sample name from the BAM header.
* You can specify manually the sample name in the `normal_sample_name` and `tumor_sample_name` fields.

Other input files can be defined for each sample, depending on the analysis starting point:
* `vcf`: the analysis will start with the annotation of vcf file with Funcotator
* `maf`: input maf files will be processed with the generation of summary plots and tables for analysis results interpretation.

An example of sample record in ``config/config.yaml`` is reported below:
```
patientA:
  normal_sample_name:
    - normalname
  tumor_sample_name:
    - tumorname
  normal_bam:
    - path/to/normal.bam
  tumor_bam:
    - path/to/tumor.bam
  vcf:
    - path/to/tumor_normal.vcf.gz
  maf:
    - path/to/sample.maf
```
Remember to delete the empty sections in a sample, i.e. if you have only a `maf` file, the sample record would look like:
```
patientA:
  maf:
    - path/to/sample.maf
```
