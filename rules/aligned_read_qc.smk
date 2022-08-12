import os

rule fastqc:
    input:
         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}.bam",

    output:
          html=config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}_fastqc.html",
          zip=config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}_fastqc.zip",
          # zip=config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}.{prefix}_fastqc.zip",
          # log:
          #      config["dir_logs"] + "aligned_reads/{cohort}/{sample}_{suffix}.GRCh38.CCS.fastqc.logs"
          # benchmark:
          #        config["dir_logs"] + "aligned_reads/{sample}_{prefix}.GRCh38.CCS.fastqc.tsv"
    threads: config["threads"]["fastqc"]["cpus"]
    params:
          extra="",
    run:
        tmp_dir = str(output.html).rstrip("_fastqc.html") + "_tmp"
        shell("mkdir -p {tmp_dir}")
        shell("{fastqc} {params.extra} -t {threads} "
              "--outdir {tmp_dir} {input.bam} "
              )


        def basename_without_ext(file_path):
            """Returns basename of file path, without the file extension."""
            base = os.path.basename(file_path)
            # split_ind = 2 if base.endswith(".bam") else 1
            base = ".".join(base.split(".")[:-1])
            return base


        output_base = basename_without_ext(input.bam)
        html_path = os.path.join(tmp_dir, output_base + "_fastqc.html")
        zip_path = os.path.join(tmp_dir, output_base + "_fastqc.zip")
        if str(output.html) != html_path:
            shell("mv {html_path} {output.html}")
        if str(output.zip) != zip_path:
            shell("mv {zip_path} {output.zip}")

rule qualimap:
    input:
         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{suffix}.bam",
    output:
          # rpt=directory(config["dir_reports"] + "aligned_reads/qualimap/{sample}/{sample}.{prefix}"),
          rpt=directory(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}", ),
    threads:
           config["threads"]["qualimap"]["cpus"]
    run:
        shell("{qualimap} --java-mem-size=100G bamqc -nt {threads} -bam {input.bam} "
              "-outdir {output.rpt} -nr 5000 -outformat PDF:HTML")
rule mosdepth:
    input:
         bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam",
    output:
          rpt=directory(config["dir_reports"] + "aligned_reads/mosdepth/{sample}/{sample}.{prefix}"),
    threads:
           config["threads"]["mosdepth"]["cpus"]
    run:
        shell("{mosdepth} -t {threads} -n {output} {input.bam}")


# rule multiqc_aligned_reads:
#     input:
#          expand(config["dir_reports"] + "aligned_reads/qualimap/{sample}/{sample}.{{prefix}}", sample=config["samples"])
#     output:
#           html=config["dir_reports"] + "aligned_reads/" + config["project"] + ".{prefix}.qualimap.multiqc.html",
#     run:
#         out_dir = "/".join(str(output.html).split("/")[:-1])
#         file_base = str(output.html).split("/")[-1].rstrip(".html")
#         info = file_base.rstrip(".qualimap.multiqc")
#         comment = "Create by Peng Jia (pengjia@stu.xjtu.edu.cn)"
#         shell("{multiqc} -o {out_dir} -n {file_base} -i '{info}' -b '{comment}' {input}")


def get_multiqc_input(wildcard):
    qualimap_dir = expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}",
                          cohort=wildcard.cohort, suffix=wildcard.suffix,
                          sample=config["samples_info"][wildcard.cohort].keys())
    fastqc_zip = expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}_fastqc.zip",
                        cohort=wildcard.cohort, suffix=wildcard.suffix,
                        sample=config["samples_info"][wildcard.cohort].keys())
    fastqc_html = expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}_fastqc.html",
                         cohort=wildcard.cohort, suffix=wildcard.suffix,
                         sample=config["samples_info"][wildcard.cohort].keys())
    return qualimap_dir + fastqc_zip + fastqc_html


rule multiqc_reads_qc:
    input:
         get_multiqc_input
         # "lslsl"
         # aligned=expand(config["dir_reports"] + "{{cohort}}/aligned_reads/{{cohort}}.{sample}.{{suffix}}.qualimap",
         #                sample=config[get_wildcard_names("cohort")]),
         # raw=expand(config["dir_reports"] + "aligned_reads/fastqc/{sample}/{sample}.{{prefix}}_fastqc.zip",
         #            sample=config["samples"])
    output:
          html=config["dir_reports"] + "final/{cohort}.{suffix}.multiqc.html",
    priority: 49

    run:
        out_dir = "/".join(str(output.html).split("/")[:-1])
        file_base = str(output.html).split("/")[-1].rstrip(".html")
        info = file_base.rstrip(".multiqc")
        comment = "Create by Peng Jia (pengjia@stu.xjtu.edu.cn)"
        shell("{multiqc} -o {out_dir} -n {file_base} -i '{info}' -b '{comment}' {input} ")


def get_multiqc_input_all(wildcard):
    qualimap_dir=[]
    fastqc_zip=[]
    fastqc_html=[]
    for cohort, cohort_info in config["samples_info"].items():
        qualimap_dir.extend(expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}",
                              cohort=cohort, suffix=wildcard.suffix,
                              sample=cohort_info.keys()))
        fastqc_zip.extend( expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}_fastqc.zip",
                            cohort=cohort, suffix=wildcard.suffix,
                            sample=cohort_info.keys()))
        fastqc_html.extend( expand(config["dir_reports"] + "{cohort}/aligned_reads/{cohort}.{sample}.{suffix}_fastqc.html",
                             cohort=cohort, suffix=wildcard.suffix,
                             sample=cohort_info.keys()))
    return qualimap_dir + fastqc_zip + fastqc_html

# def get_vcfs_for_deepvariant_merge_vcf_all(wildcards):
#     gvcf = []
#     vcf = []
#     for cohort, cohort_info in config["samples_info"].items():
#         this_gvcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
#                                                     "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.g.vcf.gz",
#                            cohort=cohort, ref_name=wildcards.ref_name, suffix=wildcards.suffix,
#                            sample=cohort_info.keys())
#         gvcf.extend(this_gvcf)
#
#         this_vcf = expand(config["dir_variants"] + "{cohort}/{sample}/deepvariant/"
#                                                    "{cohort}.{sample}.{ref_name}.{suffix}.deepvariant.raw.vcf.gz",
#                           cohort=cohort, ref_name=wildcards.ref_name, suffix=wildcards.suffix,
#                           sample=cohort_info.keys())
#         vcf.extend(this_vcf)
#     return {"gvcf_gz": gvcf, "vcf_gz": vcf}
rule multiqc_reads_qc_all:
    input:
         get_multiqc_input_all
         # "lslsl"
         # aligned=expand(config["dir_reports"] + "{{cohort}}/aligned_reads/{{cohort}}.{sample}.{{suffix}}.qualimap",
         #                sample=config[get_wildcard_names("cohort")]),
         # raw=expand(config["dir_reports"] + "aligned_reads/fastqc/{sample}/{sample}.{{prefix}}_fastqc.zip",
         #            sample=config["samples"])
    output:
          html=config["dir_reports"] + "final/"+config["project_name"]+".{suffix}.multiqc.html",
    priority: 49

    run:
        out_dir = "/".join(str(output.html).split("/")[:-1])
        file_base = str(output.html).split("/")[-1].rstrip(".html")
        info = file_base.rstrip(".multiqc")
        comment = "Create by Peng Jia (pengjia@stu.xjtu.edu.cn)"
        shell("{multiqc} -o {out_dir} -n {file_base} -i '{info}' -b '{comment}' {input} ")
