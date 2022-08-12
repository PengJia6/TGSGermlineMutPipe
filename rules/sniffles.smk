# rule samtools_chrom:
#     input:
#         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
#     output:
#         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.bam",
#     wildcard_constraints:
#             contig = "|".join(config["ref"]["GRCh38"]["avaliable"]),
#             suffix="HiFi.minimap2"
#     run:
#         shell("{samtools} view -hb {input.bam} {wildcards.contig} >{output}")
#
# rule samtools_md:
#     input:
#         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.bam",
#         ref=get_ref,
#     output:
#         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.md.bam",
#     wildcard_constraints:
#         contig="|".join(config["ref"]["GRCh38"]["avaliable"]),
#         suffix="HiFi.minimap2"
#     run:
#         shell("{samtools} calmd -b -@ 2 {input.bam} {input.ref} > {output.bam} ")
#
#
# def get_md_chrom_bam(wildcards):
#     return expand(config["dir_aligned_reads"] +
#                   "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.md.bam",
#         contig=config["ref"][wildcards.ref_name]["avaliable"],
#         cohort=wildcards.cohort,sample=wildcards.sample,ref_name=wildcards.ref_name,suffix=wildcards.suffix)
#
#
# rule samtools_merged_md:
#     input:
#        get_md_chrom_bam
#     output:
#         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.md.bam",
#     run:
#         shell("{samtools} merge -O BAM {output} {input}")


rule sniffles:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
    # sindex=config["ref"]["fasta"] + ".fai"
    output:
        vcf=config[
                "dir_variants"] + "{cohort}/{sample}/sniffles/{cohort}.{sample}.{ref_name}.{suffix}.sniffles.raw.vcf"
    params:
        extra="",
    log:
        config["dir_logs"] + "sniffles/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.sniffles.log"
    benchmark:
        config["dir_logs"] + "sniffles/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.sniffles.tsv"
    threads: config["threads"]["sniffles"]["cpus"]
    run:
        workdir = str(output.vcf)[:-4] + "_tmp"
        shell("mkdir -p {workdir}")
        shell("{sniffles} -s 3 -m {input.bam} -v {output.vcf} -t {threads} 2>{log} 1>{log}")


# shell("{bcftools} view -Oz -o {output.vcf} {vcf}")


def get_vcfs_for_sniffles_merge_vcf(wildcards):
    vcf = expand(
        config["dir_variants"] + "{cohort}/{sample}/sniffles/{cohort}.{sample}.{ref_name}.{suffix}.sniffles.raw.vcf",
        cohort=wildcards.cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,
        sample=config["samples_info"][wildcards.cohort].keys())
    return {"vcf_gz": vcf}


rule SURVIVOR_merge_vcf:
    input:
        unpack(get_vcfs_for_sniffles_merge_vcf)
    output:
        vcf=config["dir_variants"] + "{cohort}/{cohort}.{ref_name}.{suffix}.sniffles.merged.raw.vcf"
    params: ""
    priority: 50
    run:
        if len(input.vcf_gz) == 1:
            shell("cp {input.vcf_gz} {output.vcf}")
        else:
            vcf_files = str(output.vcf)[:-3] + "fofn"
            file = open(vcf_files,"w")
            for i in list(input.vcf_gz):
                file.write(f"{i}\n")
            file.close()
            shell("{SURVIVOR} merge {vcf_files} 1000 1 1 -1 -1 -1 {output.vcf}")

rule sniffles_force:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
        vcf=config["dir_variants"] + "{cohort}/{cohort}.{ref_name}.{suffix}.sniffles.merged.raw.vcf"
    # sindex=config["ref"]["fasta"] + ".fai"
    output:
        vcf=config[
                "dir_variants"] + "{cohort}/{sample}/sniffles/{cohort}.{sample}.{ref_name}.{suffix}.sniffles.raw2.vcf"
    params:
        extra="",
    log:
        config["dir_logs"] + "sniffles/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.sniffles_force.log"
    benchmark:
        config["dir_logs"] + "sniffles/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.sniffles.force.tsv"
    threads: config["threads"]["sniffles_force"]["cpus"]
    run:
        workdir = str(output.vcf)[:-4] + "_tmp"
        shell("mkdir -p {workdir}")
        shell("{sniffles} -s 3 -m {input.bam} -v {output.vcf} --Ivcf {input.vcf} -t {threads} 2>{log} 1>{log}")


def get_vcfs_for_sniffles_merge_vcf(wildcards):
    vcf = expand(
        config["dir_variants"] + "{cohort}/{sample}/sniffles/{cohort}.{sample}.{ref_name}.{suffix}.sniffles.raw2.vcf",
        cohort=wildcards.cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,
        sample=config["samples_info"][wildcards.cohort].keys())
    return {"vcf_gz": vcf}


rule SURVIVOR_merge_vcf2:
    input:
        unpack(get_vcfs_for_sniffles_merge_vcf)
    output:
        vcf=config["dir_variants"] + "{cohort}/{cohort}.{ref_name}.{suffix}.sniffles.raw.vcf"
    params: ""
    priority: 50
    run:
        if len(input.vcf_gz) == 1:
            shell("cp {input.vcf_gz} {output.vcf}")
        else:
            vcf_files = str(output.vcf)[:-3] + "fofn"
            file = open(vcf_files,"w")
            for i in list(input.vcf_gz):
                file.write(f"{i}\n")
            file.close()
            shell("{SURVIVOR} merge {vcf_files} 1000 -1 1 -1 -1 -1 {output.vcf}")


rule SURVIVOR_merge_vcf2_gzip:
    input:
        vcf=config["dir_variants"] + "{cohort}/{cohort}.{ref_name}.{suffix}.sniffles.raw.vcf"
    output:
        vcf=config["dir_variants"] + "final/{cohort}.{ref_name}.{suffix}.sniffles.raw.vcf.gz"
    run:
        shell("{bcftools} sort -Oz -o {output.vcf} {input.vcf} ")


# input_vcfs = []
# input_dir = "/".join(input.gvcf_gz[0].split("/")[:-4])
# for item in input.gvcf_gz:
#     input_vcfs.append("/".join(item.split("/")[-4:]))
# inputs_str = " ".join(["/input/" + item for item in input_vcfs])
#
# shell("docker run -v {input_dir}:/input quay.io/mlin/glnexus:v1.2.7 /usr/local/bin/glnexus_cli "
#       "--config DeepVariantWGS "
#       "{inputs_str} |"
#       "{bcftools} view -Oz -o {output.vcf}")
#
#
# def get_vcfs_for_deepvariant_merge_vcf_all(wildcards):
#     gvcf = []
#     vcf = []
#     for cohort, cohort_info in config["samples_info"].items():
#         this_gvcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
#                                                     "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
#             cohort=cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,
#             sample=cohort_info.keys())
#         gvcf.extend(this_gvcf)
#
#         this_vcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
#                                                    "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.vcf.gz",
#             cohort=cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,
#             sample=cohort_info.keys())
#         vcf.extend(this_vcf)
#     return {"gvcf_gz": gvcf, "vcf_gz": vcf}
#
#
# rule deepvariant_merge_vcf_all:
#     input:
#         # "dlld"
#         unpack(get_vcfs_for_deepvariant_merge_vcf_all)
#     # gvcf_gz=expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/deepvariant_details/"
#     #                                         "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
#     #                sample=config["samples"]),
#     # vcf_gz=expand(config["dir_variants"] + "dv/dv_details/{sample}/{sample}.{{prefix}}.dv.raw.vcf.gz",
#     #               sample=config["samples"])
#
#     output:
#         # vcf=config["dir_variants"] + "dv/" + config["project"] + ".{prefix}.dv.raw.vcf.gz"
#         vcf=config["dir_variants"] + "final/" + config["project_name"] + ".{ref_name}.{suffix}.deepvariant.raw.vcf.gz"
#     # vcf=dir_variant + "deeptrio/{trio}/{trio}.deeptrio.merged.vcf.gz",
#     # vcf_index=dir_variant + "deeptrio/{trio}/{trio}.deeptrio.merged.vcf.gz.tbi"
#     params: ""
#     priority: 50
#     # log:
#     #    config["dir_logs"] + "dv/{prefix}.combine.log"
#     #
#     # benchmark:
#     #          config["dir_logs"] + "dv/{prefix}.combine.tsv"
#
#     run:
#         if len(input.gvcf_gz) == 1:
#             shell("cp {input.vcf_gz} {output.vcf}")
#         else:
#             input_vcfs = []
#             input_dir = "/".join(input.gvcf_gz[0].split("/")[:-4])
#             for item in input.gvcf_gz:
#                 input_vcfs.append("/".join(item.split("/")[-4:]))
#             inputs_str = " ".join(["/input/" + item for item in input_vcfs])
#
#             shell("docker run -v {input_dir}:/input quay.io/mlin/glnexus:v1.2.7 /usr/local/bin/glnexus_cli "
#                   "--config DeepVariantWGS "
#                   "{inputs_str} |"
#                   "{bcftools} view -Oz -o {output.vcf}")
