# Snakemake workflow: Musta

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.15.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/musta.svg?branch=master)](https://travis-ci.org/snakemake-workflows/musta)

This workflow performs variant calling and annotation for tumoral and matched germinal samples.

## Authors

* Matteo Massidda (@massiddaMT)
* Rossano Atzeni (@ratzeni)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/).
ADD TO CATALOG
If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## TEMP INSTRUCITONS
Create a virtual environment with the command:
```commandline
mamba create -c bioconda -c conda-forge --name snakemake snakemake=6.15
```
and activate it:
```commandline
conda activate snakemake
```
You can then perform a dry-run:
```commandline
snakemake -n --use-conda -p --cores all
```
