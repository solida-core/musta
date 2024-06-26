#################
import errno
import multiprocessing
import os.path
import psutil
import yaml as yaml
from snakemake.utils import validate


report: "../report/workflow.rst"


validate(config, schema="../schemas/config.schema.yaml")

with open(config["samples"], "r") as file:
    samples_master = yaml.load(file, Loader=yaml.FullLoader)
    samples_master.keys()
    samples = []
    for sample in list(samples_master.keys()):
        if config.get("run").get("call"):
            if "tumor_bam" in samples_master[sample]:
                samples.append(sample)
            else:
                pass
        elif config.get("run").get("annotate"):
            if "vcf" in samples_master[sample]:
                samples.append(sample)
            else:
                pass
        elif config.get("run").get("analysis"):
            if "maf" in samples_master[sample]:
                samples.append(sample)
            else:
                pass
        else:
            if "tumor_bam" in samples_master[sample]:
                samples.append(sample)
            else:
                pass
d = samples_master
samples_master = {k: v for k, v in d.items() if k in samples}


def total_physical_mem_size():
    mem = psutil.virtual_memory()
    return mem.total


def cpu_count():
    return multiprocessing.cpu_count()


def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT) + " (path must be absolute)",
            filepath,
        )
    return filepath


def conservative_cpu_count(reserve_cores=1, max_cores=5):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)


def tmp_path(path=""):
    """
    if does not exists, create path and return it. If any errors, return
    default path
    :param path: path
    :return: path
    """
    default_path = os.getenv("TMPDIR", config.get("paths").get("tmp_dir"))
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                return default_path
        return path
    return default_path


def java_params(
    tmp_dir="",
    percentage_to_preserve=20,
    stock_mem=1024**3,
    stock_cpu=2,
    multiply_by=2,
):
    """
    Set Java params
    :param tmp_dir: path to tmpdir
    :param percentage_to_preserve: percentage of resources to preserve
    :param stock_mem: min memory to preserve
    :param stock_cpu: min cpu to preserve
    :param multiply_by: multiply base resource by this param
    :return: string to return to configure java environments
    """

    def bytes2human(n):
        # http://code.activestate.com/recipes/578019
        # >>> bytes2human(10000)
        # '9.8K'
        # >>> bytes2human(100001221)
        # '95.4M'
        symbols = ("K", "M", "G", "T", "P", "E", "Z", "Y")
        prefix = {}
        for i, s in enumerate(symbols):
            prefix[s] = 1 << (i + 1) * 10
        for s in reversed(symbols):
            if n >= prefix[s]:
                value = float(n) / prefix[s]
                return "%.0f%s" % (value, s)
        return "%sB" % n

    def preserve(resource, percentage, stock):
        preserved = resource - max(resource * percentage // 100, stock)
        return preserved if preserved != 0 else stock

    # def preserve(resource, percentage, stock):
    #     return resource - max(resource * percentage // 100, stock)

    params_template = "'-Xms{} -Xmx{} -XX:ParallelGCThreads={} " "-Djava.io.tmpdir={}'"

    mem_min = 1024**3 * 2  # 2GB

    mem_size = preserve(total_physical_mem_size(), percentage_to_preserve, stock_mem)
    cpu_nums = preserve(cpu_count(), percentage_to_preserve, stock_cpu)
    tmpdir = tmp_path(tmp_dir)

    return params_template.format(
        bytes2human(mem_min).lower(),
        bytes2human(max((mem_size // cpu_nums) * multiply_by, mem_min)).lower(),
        min(cpu_nums, multiply_by),
        tmpdir,
    )


def references_abs_path(ref="ref"):
    references = config.get(ref)
    basepath = expand_filepath(references["basepath"])
    provider = references["provider"]
    release = references["release"]

    return [os.path.join(basepath, provider, release)]


def resolve_single_filepath(basepath, filename):
    return os.path.join(basepath, filename)


def get_name(inputfile):
    with open(inputfile, "r") as file:
        data = file.read().rstrip()
    return data


def get_tumoral_bam(wildcards):
    with open(config["samples"], "r") as file:
        samples_master = yaml.load(file, Loader=yaml.FullLoader)
        samples_master.keys()
    return samples_master[wildcards.sample]["tumor_bam"][0]


def get_normal_bam(wildcards):
    with open(config["samples"], "r") as file:
        samples_master = yaml.load(file, Loader=yaml.FullLoader)
        samples_master.keys()
    # print(wildcards.sample)
    return samples_master[wildcards.sample]["normal_bam"][0]

def get_vep_genome_version(version):
    if version in ['hg19', 'hg38']:
        return 'GRCh37' if version in 'hg19' else 'GRCh38'
    return version

def select_filtered(wildcards):
    with open(config["samples"], "r") as file:
        samples_master = yaml.load(file, Loader=yaml.FullLoader)
    if not samples_master[wildcards.sample]["normal_bam"]:
        return rules.filter_mutect_tumoronly.output.vcf
    else:
        return rules.filter_mutect.output.vcf


## functions for pipeline starting from vcf


def get_annotation_input(wildcards):
    with open(config["samples"],"r") as file:
        samples_master = yaml.load(file,Loader=yaml.FullLoader)
    if not samples_master[wildcards.sample]["vcf"]:
        return rules.somaticseq_hold_on.output.snvs
    else:
        return samples_master[wildcards.sample]["vcf"][0]

#
# def get_vcf_list(wildcards):
#     with open(config["samples"], "r") as file:
#         samples_master = yaml.load(file, Loader=yaml.FullLoader)
#         samples_master.keys()
#
#     return samples_master[wildcards.sample]["vcf"][0] if "vcf" in samples_master[wildcards.sample] else None


def get_maf_file_input():
    if config.get("run").get("interpret"):
        return lambda wildcards: get_maf_list(wildcards)
    else:
        return expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "results/annotation/funcotator/{sample}_funcotated.maf",
            ),
            sample=list(samples_master.keys()),
        )


def get_maf_list(wildcards):
    with open(config["samples"], "r") as file:
        samples_master = yaml.load(file, Loader=yaml.FullLoader)
    return (p[1].get("maf")[0] for p in samples_master.items() if "maf" in p[1])


def get_normalname(wildcards):
    with open(config["samples"], "r") as file:
        samples_master = yaml.load(file, Loader=yaml.FullLoader)
        samples_master.keys()
    return samples_master[wildcards.sample]["normal_sample_name"][0]


def get_tumorname(wildcards):
    with open(config["samples"], "r") as file:
        samples_master = yaml.load(file, Loader=yaml.FullLoader)
        samples_master.keys()
    return samples_master[wildcards.sample]["tumor_sample_name"][0]


def ensure_dir(path, force=False):
    try:
        if force and os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def exist_dir(path, delete=False):
    try:
        if delete and os.path.exists(path):
            shutil.rmtree(path)
    # os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def resolve_results_filepath(basepath, outname):
    return os.path.join(basepath, outname)

def generate_flag(flag_name, input_file):
    if input_file:
        return "--{} {}".format(flag_name, input_file)
    else:
        return ""
