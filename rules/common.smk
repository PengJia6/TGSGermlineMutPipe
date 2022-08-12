rule bam2fastqgz:
    input:
         "{prefix}.bam"
    output:
          "{prefix}.fq.gz"
    threads: config["threads"]["bam2fastqgz"]["cpus"]
    run:
        shell("{samtools} fastq -0 {output} -@ {threads} {input}")

rule bam2fastqgz2:
    input:
         "{prefix}.bam"
    output:
          "{prefix}.fastq.gz"
    threads: config["threads"]["bam2fastqgz2"]["cpus"]
    run:
        shell("{samtools} fastq -0 {output} -@ {threads} {input}")
#
# rule sample_bam_merge:
#     input:
#          get_ngmlr_merge_bams
#     output:
#           bam=dir_aligned_reads + "{sample}/{sample}.GRCh38.CCS.{aligner}.addrg.bam",
#     log:
#        dir_logs + "{sample}.GRCh38.CCS.{aligner}.logs"
#     benchmark:
#              dir_logs + "{sample}.GRCh38.CCS.{aligner}.tsv"
#
#     threads: config["threads"]["sample"]
#     run:
#         if len(input) < 2:
#             shell("ln -sr {input} {output.bam} ")
#             shell("sleep 1")
#             shell("touch -h {output.bam} 1>>{log} 2>>{log}")
#             shell("echo only one bam, make soft link 1>>{log} 2>>{log}")
#         else:
#             shell("{samtools} merge -@ {threads} {output.bam} {input} 2>{log} 1>{log}")

rule bam_index:
    input:
         "{prefix}.bam"
    output:
          "{prefix}.bam.bai"
    threads: config["threads"]["bam_index"]["cpus"]
    shell:
         "{samtools} index  -@ {threads} {input}"

rule tabix:
    input:
         "{prefix}.vcf.gz"
    output:
          "{prefix}.vcf.gz.tbi"
    run:
        shell("{tabix} -f {input}")
