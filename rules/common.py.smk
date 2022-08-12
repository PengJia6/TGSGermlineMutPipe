import yaml


def get_more_samples_info(conf_samples):
    file_samples = open(conf_samples, "r")
    samples_info = yaml.safe_load(file_samples)
    # print(samples_info)
    for cohort, cohort_info in samples_info.items():
        for sample_name, info in cohort_info.items():
            if "path" not in info:
                print("E: There is not 'path' info in sample configure file!")
                return {}
            this_sub_samples = {}
            for one_path in info["path"]:
                if not (one_path.endswith("fastq.gz") or one_path.endswith("fq.gz")):
                    print("E: The path {} not end with fastq.gz or fq.gz!".format(one_path))
                    return {}
                if len(one_path) > 2:
                    sub_name = one_path.split("/")[-1].rstrip("fq.gz").rstrip("fastq.gz")
                    this_sub_samples[sub_name] = one_path
            if len(this_sub_samples) < 1:
                print("E: There is not available path in sample {}".format(sample_name))
                return {}
            samples_info[cohort][sample_name]["subsamples"] = this_sub_samples
    return samples_info


# def get_one_sample_fq(wildcards):
#     return samples_info[wildcards.sample]["subsamples"][wildcards.subsample]


def add_cluster(conf_cluster):
    file_config_cluster = open(conf_cluster, "r")
    conf_cluster = yaml.safe_load(file_config_cluster)
    return conf_cluster


def get_ref(wildcards):
    return config["ref"][wildcards.ref_name]["fasta"]


#
# def get_one_unmaped_fq(wildcards):
#     sample_name = str(wildcards.sample)
#     mysubsample = str(wildcards.subsample)
#     raw_bam_list = open(dir_raw_reads + "{sample}/{sample}.CCS.fofn".format(sample=sample_name), "r").readlines()
#     this_samples = {}
#     for line in raw_bam_list:
#         if len(line) < 3: continue
#         subsample = line[:-1].split("/")[-1].rstrip("fastq.gz").rstrip("fq.gz")
#         this_samples[subsample] = line.rstrip("\n")
#     if mysubsample in this_samples:
#         return this_samples[mysubsample]
#     else:
#         return "error input"
samples_info = get_more_samples_info(config["path_samples"])  # see "rules/common.smk"
config["threads"] = add_cluster(config["path_cluster_config"])  # see "rules/common.smk"
# config["contigs"] = config["ref"]["autosomes"] + config["ref"]["sex_chrom"] + config["ref"]["mit_chrom"]
config["cohorts"] = list(samples_info.keys())
config["samples"] = list(samples_info.keys())
config["samples_info"] = samples_info
