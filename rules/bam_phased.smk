# rule read_phasing:
#     input:
#          bam=config["dir_tmp"] + "align/{sample}/{sample}.GRCh38.CCS.{aligner}.bam",
#          bai=config["dir_tmp"] + "align/{sample}/{sample}.GRCh38.CCS.{aligner}.bam.bai",
#          ref=config["ref"]["fasta"],
#          vcf_gz=config["dir_variants"] + "dv/{sample}/{sample}.GRCh38.CCS.{aligner}.dv.phased.vcf.gz",
#          vcf_gz_index=config["dir_variants"] + "dv/{sample}/{sample}.GRCh38.CCS.{aligner}.dv.phased.vcf.gz.tbi",
#
#     output:
#           bam=config['dir_aligned_reads'] + "{sample}/{sample}.GRCh38.CCS.{aligner}.dv_phased.bam",
#     threads: config["threads"]["read_phasing"]
#     log:
#        config["dir_logs"] + "read_phasing/{sample}.GRCh38.CCS.{aligner}.dv.read_phased.logs"
#     benchmark:
#              config["dir_logs"] + "read_phasin/{sample}.GRCh38.CCS.{aligner}.dv.read_phased.tsv"
#     run:
#         shell("{whatshap} haplotag --output {output.bam} --reference {input.ref} "
#               "{input.vcf_gz} {input.bam} 2>{log} 1>{log}")

rule read_phasing:
    input:
         bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam",
         bai=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam.bai",
         ref=config["ref"]["fasta"],
         vcf_gz=config["dir_variants"] +
                "{caller_p}/{caller_p}_details/{sample}/{sample}.{prefix}.{caller_p}.raw.phased4read.vcf.gz",
         vcf_gz_index=config["dir_variants"] +
                      "{caller_p}/{caller_p}_details/{sample}/{sample}.{prefix}.{caller_p}.raw.phased4read.vcf.gz.tbi"

    output:
          bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.{caller_p}_phased.bam",
    threads: config["threads"]["read_phasing"]["cpus"]
    log:
       config["dir_logs"] + "read_phasing/{sample}/{sample}.{prefix}.{caller_p}.log"
    benchmark:
             config["dir_logs"] + "read_phasing/{sample}/{sample}.{prefix}.{caller_p}.tsv"
             # config["dir_logs"] + "read_phasin/{sample}.GRCh38.CCS.{aligner}.dv.read_phased.tsv"
    run:
        shell("{whatshap} haplotag --output {output.bam} --reference {input.ref} "
              "{input.vcf_gz} {input.bam} 2>{log} 1>{log}")
