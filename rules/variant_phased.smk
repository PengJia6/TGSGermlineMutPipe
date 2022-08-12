def get_bams(wildcards):
    return config["bams_hq"][wildcards.aligner]


ruleorder: var_phased > var_phased_for_read_phased

rule var_phased:
    input:
         vcf="{prefix_all}.{aligner}.{suffix}.vcf.gz",
         bams=get_bams,
         ref=config["ref"]["fasta"]
    output:
          "{prefix_all}.{aligner}.{suffix}.phased.vcf.gz",
          # log: ""
          # benchmark: ""
    threads: config["threads"]["var_phased"]["cpus"]
    run:
        shell("{whatshap} phase --output {output} --reference {input.ref} {input.vcf} "
              "{input.bams} ")

##TODO
rule var_phased_for_read_phased:
    input:
         bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam",
         bai=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam.bai",
         ref=config["ref"]["fasta"],
         vcf_gz=config[
                    "dir_variants"] + "{caller_p}/{caller_p}_details/{sample}/{sample}.{prefix}.{caller_p}.raw.vcf.gz",
         vcf_gz_index=config[
                          "dir_variants"] + "{caller_p}/{caller_p}_details/{sample}/{sample}.{prefix}.{caller_p}.raw.vcf.gz.tbi",
    output:
          # =config["dir_variants"] + "{prefix}.dv.phased.vcf.gz",
          vcf_gz=config[
                     "dir_variants"] + "{caller_p}/{caller_p}_details/{sample}/{sample}.{prefix}.{caller_p}.raw.phased4read.vcf.gz"

    log:
       config["dir_logs"] + "{caller_p}_phased/{sample}.{prefix}.{caller_p}.phased.logs"
    benchmark:
             config["dir_logs"] + "{caller_p}_phased/{sample}.{prefix}.{caller_p}.phased.tsv"
    threads: config["threads"]["var_phased_for_read_phased"]["cpus"]
    run:
        shell("{whatshap} phase --output {output.vcf_gz} --reference {input.ref} {input.vcf_gz} "
              "{input.bam} 2>{log} 1>{log}")
