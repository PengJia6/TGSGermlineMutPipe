import re
from pathlib import Path
import yaml

# shell.prefix("set -o pipefail; umask 002; ")  # set g+w
configfile: "conf/reference.yaml"  # reference information
configfile: "conf/config.yaml"  # reference information
# configfile: "conf/samples_tmp2.yaml"  # reference information
config["path_cluster_config"] = "conf/cluster.yaml"
config["path_samples"] = "conf/samples.yaml"
dir_data = config["dir_data"]
config["dir_raw_data"] = dir_data + "rawdata/"
config["dir_aligned_reads"] = dir_data + "aligned_reads/"
config["dir_reports"] = dir_data + "reports/"
config["dir_tmp"] = dir_data + "tmp/"
config["dir_variants"] = dir_data + "variants/"
config["dir_logs"] = dir_data + "logs/"
# conf_samples = "conf/samples_tmp2.yaml"
# refs = config["ref"]
# aligners = config["aligners"]
# callers = config["callers"]

### dir for this project

# path_ref = "/home/DATA/Chinese_Quartet/ref_based_analysis/data_pre/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
include: "rules/common.py.smk"
include: "rules/software.smk"

# samples = ["LCL5", "LCL6", "LCL7", "LCL8"]
# for sample in samples:
#     info = open(dir_data + "rawdata/HiFi/fofn/{}.HiFi.fofn".format(sample), "r").readlines()
#     # thisinfo = {i.rstrip(".bam\n").rstrip(".fq.gz").rstrip(".fastq.gz").split("/")[-1]: i.rstrip("\n") for i in info}
#     samples_info[sample] = [i.rstrip("\n") for i in info]
samples_info = config["samples_info"]
# print(samples_info)
# print([item for cohort, info in config["samples_info"].items() for item in info])
wildcard_constraints:
    caller="|".join(config["callers"] + config["svcallers"]),
    aligner="|".join(config["aligners"]),
    # var_type="SNV|INDEL|SV|CNV",
    cohort="|".join((samples_info.keys())),
    ref_name="|".join(config["refs"]),
    sample="|".join([item for cohort, info in config["samples_info"].items() for item in info]),
# caller_p="|".join(confi
bam_prefix = "{ref_name}.HiFi.{aligner}"
if config["bam_left_align"]:
    bam_prefix = bam_prefix + "." + "leftAlign"
if config["bam_phasing"]:
    bam_phasing_caller = config["bam_phasing_caller"]
    bam_prefix = bam_prefix + "." + bam_phasing_caller + "_phased"

targets = []
# bams_hq = {}
# bams_hq_dict = {}
for ref_name in config["refs"]:
    for aligner in config["aligners"]:
        # bams_hq_dict[aligner] = {}
        # bams_hq_this_aligner = []
        this_prefix = bam_prefix.format(aligner=aligner,ref_name=ref_name)
        # targets.append(
        #     config["dir_reports"] + config["project_name"] + ".{prefix}.multiqc.html".format(
        #         prefix=this_prefix,
        #     ))
        # if "deepvariant" in config["callers"]:
        #     targets.append((config["dir_variants"] +
        #                     "final/{cohort}.{prefix}.{caller}.raw.vcf.gz".format(
        #                         prefix=this_prefix, cohort=config["project_name"], caller="deepvariant")))
        for cohort, cohort_info in samples_info.items():
            # for sample in cohort_info:
            # hq_bam = config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{prefix}.phased.bam.bai".format(
            #     cohort=cohort,
            #     sample=sample,
            #     prefix=this_prefix)
            # targets.append(hq_bam)
            # bams_hq_this_aligner.append(hq_bam)
            # bams_hq_dict[aligner][sample] = hq_bam
            for caller in config["svcallers"] + config["callers"]:
                if caller in ["sniffles", "deepvariant"]:
                    targets.append((config["dir_variants"] +
                                    "final/{cohort}.{prefix}.{caller}.raw.vcf.gz".format(
                                        prefix=this_prefix,cohort=cohort,caller=caller)))
                elif caller in ["cuteSV", "pbsv", "SVision"]:
                    for sample in cohort_info:
                        targets.append(config["dir_variants"] +
                                       "{cohort}/{sample}/{caller}/{cohort}.{sample}.{prefix}.{caller}.raw.vcf.gz".format(
                                           sample=sample,prefix=this_prefix,cohort=cohort,caller=caller))
# for sample in cohort_info:
#     targets.append(config["dir_variants"] +
#                    "{cohort}/{sample}/{caller}/{cohort}.{sample}.{prefix}.{caller}.raw.vcf.gz".format(
#                        sample=sample, prefix=this_prefix, cohort=cohort, caller=caller))
# if caller == "deepvariant":
#     targets.append((config["dir_variants"] +
#                     "final/{cohort}.{prefix}.{caller}.raw.vcf.gz".format(
#                         prefix=this_prefix, cohort=cohort, caller=caller)))
# bam_report = config["dir_reports"] + "final/{cohort}.{prefix}.multiqc.html".format(prefix=this_prefix,
#                                                                                    cohort=cohort)
# targets.append(bam_report)

# bams_hq[aligner] = bams_hq_this_aligner

# for caller in callers:
#     # if  caller == config["bam_phasing_caller"] : continue
#     # callers_prefix = bam_prefix + ".{caller}.final"
#     targets.append(config["dir_variants"] + "{}/".format(caller) + config["project"] +
#                    ".{prefix}.{caller}.marked.norm.phased.vcf.gz".format(prefix=this_prefix, caller=caller))
#
# targets.append(config["dir_reports"] + "aligned_reads/" + config["project"] +
#                ".{prefix}.qualimap.multiqc.html".format(prefix=this_prefix))
# targets.append(config["dir_reports"] + config["project"] +
#                ".{prefix}.multiqc.html".format(prefix=this_prefix))
# targets.append(config["dir_variants"] + "{}/".format(caller) + config["project"] +
#                ".{prefix}.{caller}.marked.norm.phased.vcf.gz".format(prefix=this_prefix, caller=caller))

# config["bams_hq"] = bams_hq
# config["bams_hq_dict"] = bams_hq_dict
# config["bams_prefix"] = bams_prefix
# targets.append(config["dir_reports"] + "raw_reads/" + config["project"] + ".fastqc.multiqc.html")

include: "rules/common.smk"
include: "rules/minimap2.smk"
#
# include: "rules/pbmm2.smk"
# include: "rules/bam_merge.smk"
# include: "rules/leftalign.smk"
include: "rules/deepvariant.smk"
# include: "rules/bam_phased.smk"
# include: "rules/gatk.smk"
# include: "rules/varscan.smk"
# include: "rules/vcf_normalization.smk"
# include: "rules/variant_phased.smk"
# include: "rules/raw_read_qc.smk"
include: "rules/aligned_read_qc.smk"
include: "rules/cuteSV.smk"
include: "rules/chineseQuartetphased.smk"
include: "rules/pbsv.smk"
include: "rules/sniffles.smk"
include: "rules/SVision.smk"

# targets.append(dir_aligned_reads_tmp + "")
# include: 'rules/sample_common.smk'
# # call structural variants with pbsv
# include: 'rules/sample_pbsv.smk'
# if 'pbsv_vcf' in config['sample_targets']:
#     # pbsv VCFs
#     targets.extend([f"samples/{sample}/pbsv/{sample}.{ref}.pbsv.{suffix}"
#                     for suffix in ['vcf.gz', 'vcf.gz.tbi']])
#
# # call small variants with DeepVariant
# include: 'rules/sample_deepvariant.smk'
# if 'deepvariant' in config['sample_targets']:
#     # deepvariant VCFs, gVCFs, reports, and stats
#     targets.extend([f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.{suffix}"
#                     for suffix in ['vcf.gz', 'vcf.gz.tbi', 'g.vcf.gz', 'g.vcf.gz.tbi',
#                                    'visual_report.html', 'vcf.stats.txt']])
#
# # phase small variants with WhatsHap and haplotag BAM
# include: 'rules/sample_whatshap.smk'
# if 'whatshap' in config['sample_targets']:
#     # phased VCFs, stats, phase block GTFs, and haplotagged BAMs
#     targets.extend([f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.{suffix}"
#                     for suffix in ['phased.vcf.gz', 'phased.vcf.gz.tbi', 'phased.gtf',
#                                    'phased.tsv', 'phased.blocklist',
#                                    'haplotagged.bam', 'haplotagged.bam.bai']])
#
# # genotype STRs
# include: 'rules/sample_tandem_genotypes.smk'
# if 'tandem-genotypes' in config['sample_targets']:
#     # tandem-genotypes tabular output and plots
#     targets.extend([f"samples/{sample}/tandem-genotypes/{sample}.tandem-genotypes.{suffix}"
#                     for suffix in ['txt', 'pdf']])
#
# # calculate coverage of haplotagged sample aBAM with mosdepth
# include: 'rules/sample_mosdepth.smk'
# include: 'rules/sample_gc_coverage.smk'
# if 'coverage' in config['sample_targets']:
#     # coverage from merged haplotagged aBAM
#     targets.extend([f"samples/{sample}/mosdepth/{sample}.{ref}.deepvariant.haplotagged.{suffix}"
#                     for suffix in ['mosdepth.global.dist.txt', 'mosdepth.region.dist.txt',
#                                    'mosdepth.summary.txt', 'regions.bed.gz']])
#     targets.extend([f"samples/{sample}/mosdepth/{sample}.{ref}.gc_coverage.summary.txt"])

localrules: all

rule all:
    input:
        targets
    # "aa",
    # targets + [f"{x}.md5" for x in targets]
    #
    # rule md5sum:
    #     input: "{prefix}"
    #     output: "{prefix}.md5"
    #     message: "Creating md5 checksum for {input}."

    #     shell: "md5sum {input} > {output}"
