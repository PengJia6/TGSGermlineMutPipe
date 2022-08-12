rule deepvariant:
    input:
         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
         bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
         ref=get_ref
    output:
          vcf_gz=config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
                                          "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.vcf.gz",
          gvcf_gz=config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
                                           "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
          # gvcf_gz=config["dir_variants"] + "dv/dv_details/{sample}/{sample}.{prefix}.dv.raw.g.vcf.gz"
    log:
       config["dir_logs"] + "deepvariant/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.deepvariant.logs"
       # benchmark:
       #          config["dir_logs"] + "dv/{sample}/{sample}.{prefix}.dv.tsv"
    benchmark:
             config["dir_logs"] + "deepvariant/{cohort}/{sample}/{sample}.{ref_name}.{suffix}.deepvariant.logs"

    threads: config["threads"]["deepvariant"]["cpus"]
    run:
        dir_tmp = str(output.vcf_gz).rstrip(".vcf.gz") + "_tmp"
        file_tmp = dir_tmp.split("/")[-1]
        shell("mkdir -p " + dir_tmp)
        bam_dir = "/".join(str(input.bam).split("/")[:-1])
        bam_file = str(input.bam).split("/")[-1]
        ref_dir = "/".join(str(input.ref).split("/")[:-1])
        ref_file = str(input.ref).split("/")[-1]
        output_dir = "/".join(str(output.vcf_gz).split("/")[:-1])
        output_file = str(output.vcf_gz).split("/")[-1].rstrip(".vcf.gz")

        shell('docker run '
              '-v "{bam_dir}":"/input" '
              '-v "{ref_dir}":"/ref" '
              '-v "{output_dir}":"/output" '
              'google/deepvariant:1.1.0 /opt/deepvariant/bin/run_deepvariant '
              '--model_type=PACBIO '
              '--ref=/ref/{ref_file} '
              '--reads=/input/{bam_file} '
              '--output_vcf=/output/{output_file}.vcf '
              '--output_gvcf=/output/{output_file}.g.vcf '
              '--num_shards={threads} '
              '--make_examples_extra_args min_mapping_quality=1,keep_supplementary_alignments=true '
              '--intermediate_results_dir /output/{file_tmp} 1>{log} 2>{log}')
        shell("{bcftools} view -Oz -o {output.vcf_gz} {output_dir}/{output_file}.vcf")
        shell("{bcftools} view -Oz -o {output.gvcf_gz} {output_dir}/{output_file}.g.vcf")


# shell("touch {output.vcf_gz}")

def get_vcfs_for_deepvariant_merge_vcf(wildcards):
    gvcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
                                           "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
                  cohort=wildcards.cohort, ref_name=wildcards.ref_name, suffix=wildcards.suffix,
                  sample=config["samples_info"][wildcards.cohort].keys())
    vcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
                                          "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.vcf.gz",
                 cohort=wildcards.cohort, ref_name=wildcards.ref_name, suffix=wildcards.suffix,
                 sample=config["samples_info"][wildcards.cohort].keys())
    return {"gvcf_gz": gvcf, "vcf_gz": vcf}


rule deepvariant_merge_vcf:
    input:
         # "dlld"
         unpack(get_vcfs_for_deepvariant_merge_vcf)
         # gvcf_gz=expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/deepvariant_details/"
         #                                         "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
         #                sample=config["samples"]),
         # vcf_gz=expand(config["dir_variants"] + "dv/dv_details/{sample}/{sample}.{{prefix}}.dv.raw.vcf.gz",
         #               sample=config["samples"])

    output:
          # vcf=config["dir_variants"] + "dv/" + config["project"] + ".{prefix}.dv.raw.vcf.gz"
          vcf=config["dir_variants"] + "final/{cohort}.{ref_name}.{suffix}.deepvariant.raw.vcf.gz"
          # vcf=dir_variant + "deeptrio/{trio}/{trio}.deeptrio.merged.vcf.gz",
          # vcf_index=dir_variant + "deeptrio/{trio}/{trio}.deeptrio.merged.vcf.gz.tbi"
    params: ""
    priority: 50
    # log:
    #    config["dir_logs"] + "dv/{prefix}.combine.log"
    #
    # benchmark:
    #          config["dir_logs"] + "dv/{prefix}.combine.tsv"

    run:
        if len(input.gvcf_gz) == 1:
            shell("cp {input.vcf_gz} {output.vcf}")
        else:
            input_vcfs = []
            input_dir = "/".join(input.gvcf_gz[0].split("/")[:-4])
            for item in input.gvcf_gz:
                input_vcfs.append("/".join(item.split("/")[-4:]))
            inputs_str = " ".join(["/input/" + item for item in input_vcfs])

            shell("docker run -v {input_dir}:/input quay.io/mlin/glnexus:v1.2.7 /usr/local/bin/glnexus_cli "
                  "--config DeepVariantWGS "
                  "{inputs_str} |"
                  "{bcftools} view -Oz -o {output.vcf}")


def get_vcfs_for_deepvariant_merge_vcf_all(wildcards):
    gvcf = []
    vcf = []
    for cohort, cohort_info in config["samples_info"].items():
        this_gvcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
                                                    "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
                           cohort=cohort, ref_name=wildcards.ref_name, suffix=wildcards.suffix,
                           sample=cohort_info.keys())
        gvcf.extend(this_gvcf)

        this_vcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
                                                   "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.vcf.gz",
                          cohort=cohort, ref_name=wildcards.ref_name, suffix=wildcards.suffix,
                          sample=cohort_info.keys())
        vcf.extend(this_vcf)
    return {"gvcf_gz": gvcf, "vcf_gz": vcf}


rule deepvariant_merge_vcf_all:
    input:
         # "dlld"
         unpack(get_vcfs_for_deepvariant_merge_vcf_all)
         # gvcf_gz=expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/deepvariant_details/"
         #                                         "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
         #                sample=config["samples"]),
         # vcf_gz=expand(config["dir_variants"] + "dv/dv_details/{sample}/{sample}.{{prefix}}.dv.raw.vcf.gz",
         #               sample=config["samples"])

    output:
          # vcf=config["dir_variants"] + "dv/" + config["project"] + ".{prefix}.dv.raw.vcf.gz"
          vcf=config["dir_variants"] + "final/" + config["project_name"] + ".{ref_name}.{suffix}.deepvariant.raw.vcf.gz"
          # vcf=dir_variant + "deeptrio/{trio}/{trio}.deeptrio.merged.vcf.gz",
          # vcf_index=dir_variant + "deeptrio/{trio}/{trio}.deeptrio.merged.vcf.gz.tbi"
    params: ""
    priority: 50
    # log:
    #    config["dir_logs"] + "dv/{prefix}.combine.log"
    #
    # benchmark:
    #          config["dir_logs"] + "dv/{prefix}.combine.tsv"

    run:
        if len(input.gvcf_gz) == 1:
            shell("cp {input.vcf_gz} {output.vcf}")
        else:
            input_vcfs = []
            input_dir = "/".join(input.gvcf_gz[0].split("/")[:-4])
            for item in input.gvcf_gz:
                input_vcfs.append("/".join(item.split("/")[-4:]))
            inputs_str = " ".join(["/input/" + item for item in input_vcfs])

            shell("docker run -v {input_dir}:/input quay.io/mlin/glnexus:v1.2.7 /usr/local/bin/glnexus_cli "
                  "--config DeepVariantWGS "
                  "{inputs_str} |"
                  "{bcftools} view -Oz -o {output.vcf}")

# shell("{tabix} {output.vcf}")
#
# rule deepvariant_combile_vcf:
#     input:
#          vcfs=expand(config["dir_variants"] + "dv/dv_details/{sample}/{sample}.{{prefix}}.dv.raw.vcf.gz",
#                      sample=config["samples"]),
#          vcfs_index=expand(config["dir_variants"] + "dv/dv_details/{sample}/{sample}.{{prefix}}.dv.raw.vcf.gz.tbi",
#                            sample=config["samples"]),
#     output:
#           config["dir_variants"] + "dv/" + config["project"] + ".{prefix}.dv.raw.vcf.gz"
#     threads: config["threads"]["deepvariant_combile_vcf"]["cpus"]
#     log:
#        config["dir_logs"] + "dv/{prefix}.combine.log"
#     benchmark:
#              config["dir_logs"] + "dv/{prefix}.combine.tsv"
#     run:
#         shell("{bcftools} merge -Oz -0 --threads {threads} -o {output} {input.vcfs}")
# rule deepvariant_marker_low_quailty:
#     input:
#          vcf=config["dir_variants"] + "dv/" + config["project"] + ".{prefix}.dv.raw.vcf.gz",
#          vcf_index=config["dir_variants"] + "dv/" + config["project"] + ".{prefix}.dv.raw.vcf.gz.tbi"
#     output:
#           vcf=config["dir_variants"] + "dv/" + config["project"] + ".{prefix}.dv.marked.vcf.gz"
#     threads: config["threads"]["deepvariant_marker_low_quailty"]["cpus"]
#     log:
#        config["dir_logs"] + "dv/" + config["project"] + ".{prefix}.dv_marked_low_quality.log"
#     benchmark:
#              config["dir_logs"] + "dv/" + config["project"] + ".{prefix}.dv_marked_low_quality.tsv"
#     run:
#         shell("{bcftools} filter -m + -sLowQual_PJ "
#               "-i '%QUAL>=10 && DP<300 && DP>=10 ' "
#               "-Oz -o {output.vcf} {input.vcf} ")
