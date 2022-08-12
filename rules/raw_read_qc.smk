

rule multiqc_raw_reads:
    input:
         expand(config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}.{{prefix}}_fastqc.zip",
                sample=config["samples"])
    output:
          html=config["dir_reports"] + "raw_reads/" + config["project"] + ".{prefix}.fastqc.multiqc.html",
    run:
        out_dir = "/".join(str(output.html).split("/")[:-1])
        file_base = str(output.html).split("/")[-1].rstrip(".html")
        info = file_base.rstrip(".fastqc.multiqc")
        comment = "Create by Peng Jia (pengjia@stu.xjtu.edu.cn)"
        shell("{multiqc} -o {out_dir} -n {file_base} -i '{info}' -b '{comment}' {input}")

#
# rule fastqc:
#     input:
#          fq=get_one_sample_fq,  #rules.common
#     output:
#           html=config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}_{subsample}_fastqc.html",
#           zip=config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}_{subsample}_fastqc.zip",
#     log:
#        config["dir_logs"] + "raw_reads/{sample}_{subsample}.GRCh38.CCS.fastqc.logs"
#     benchmark:
#              config["dir_logs"] + "align/{sample}_{subsample}.GRCh38.CCS.fastqc.tsv"
#     threads: config["threads"]["fastqc"]["cpus"]
#     params:
#           extra="",
#     run:
#         tmp_dir = str(output.html).rstrip(".fastqc.html") + "_tmp"
#         shell("mkdir -p {tmp_dir}")
#
#         shell("{fastqc} {params.extra} -t {threads} "
#               "--outdir {tmp_dir} {input.fq} "
#               " 2>>{log} 1>>{log}")
#
#
#         def basename_without_ext(file_path):
#             """Returns basename of file path, without the file extension."""
#             base = os.path.basename(file_path)
#             split_ind = 2 if base.endswith(".gz") else 1
#             base = ".".join(base.split(".")[:-split_ind])
#             return base
#
#
#         output_base = basename_without_ext(input.fq)
#         html_path = os.path.join(tmp_dir, output_base + "_fastqc.html")
#         zip_path = os.path.join(tmp_dir, output_base + "_fastqc.zip")
#         if str(output.html) != html_path:
#             shell("mv {html_path} {output.html}")
#         if str(output.zip) != zip_path:
#             shell("mv {zip_path} {output.zip}")
#
#
# def get_fastqc_list(wildcards):
#     subsamples = list(samples_info[wildcards.sample]["subsamples"].keys())
#     return expand(config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}_{subsample}_fastqc.zip",
#                   sample=wildcards.sample, subsample=subsamples)
#
#
# rule fastqc_report_list:
#     input:
#          get_fastqc_list
#     output:
#           config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}.fastqc.list",
#
#     run:
#         file = open(str(output), "w")
#         for i in input:
#             file.write(str(i) + "\n")
#         file.close()
#
# rule multiqc_raw_reads:
#     input:
#          expand(config["dir_reports"] + "raw_reads/fastqc/{sample}/{sample}.fastqc.list", sample=config["samples"])
#     output:
#           fastqc_list=config["dir_reports"] + "raw_reads/" + config["project"] + ".fastqc.list",
#           html=config["dir_reports"] + "raw_reads/" + config["project"] + ".fastqc.multiqc.html",
#     run:
#         shell("cat {input} > {output.fastqc_list}")
#         if os.path.exists(output.html):
#             os.system("rm {}".format(output.html))
#         out_dir = "/".join(str(output.html).split("/")[:-1])
#         file_base = str(output.html).split("/")[-1].rstrip(".html")
#         info = file_base.rstrip(".fastqc.multiqc")
#         comment = "Create by Peng Jia (pengjia@stu.xjtu.edu.cn)"
#         shell("{multiqc} -o {our_dir} -n {file_base} -i '{info}' -b '{comment}' -l {output.fastqc_list}")
