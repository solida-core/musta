## scrivi funzione per decidere l'input

def get_annotation_input():
    if config.get("run").get("annotate"):
        return(lambda wildcards: get_vcf_list(wildcards))
    else:
        return "results/{sample}_somatic_filtered_selected.vcf.gz"

def get_vcf_list(wildcards):
    with open(config["samples"],'r') as file:
        samples_master = yaml.load(file,Loader=yaml.FullLoader)
        samples_master.keys()
    # print(wildcards.sample)
    return samples_master[wildcards.sample]["vcf"][0]



rule Funcotator:
    input:
        vcf=get_annotation_input()
    output:
        vcf=report(
            "results/annotation/funcotator/{sample}_funcotated.maf",
            caption="../report/vcf.rst",
            category="Annotation",
        )
    params:
        custom=java_params(tmp_dir=config.get("processing").get("tmp_dir"),multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path("ref"),config.get("ref").get("fasta")),
        intervals=config.get("processing").get("interval_list"),
        resources=config.get("params").get("gatk").get("Funcotator")
    log:
        "logs/gatk/Funcotator/{sample}.funcotator.log"
    conda:
       "../envs/gatk.yaml"
    threads:
        conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir = config.get("processing").get("tmp_dir")
    shell:
        "gatk Funcotator "
        "--java-options {params.custom} "
        "-R {params.genome} "
        "-L {params.intervals} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "--data-sources-path {params.resources} "
        "--output-file-format MAF "
        "--ref-version hg19 "
        ">& {log} "


