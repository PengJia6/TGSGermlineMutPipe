rule samtools_mpileup:
    input:
         bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam",
         bai=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam.bai",
         ref=config["ref"]["fasta"],
         sindex=config["ref"]["fasta"] + ".fai"
    output:
          config["dir_variants"] + "varscan/varscan_details/{sample}/{sample}.{prefix}.mpileup"
    params:
          extra="",
          dp=5
    log:
       config["dir_logs"] + "varscan/{sample}/{sample}.{prefix}.samtools_mpileup.call.log"
    benchmark:
             config["dir_logs"] + "varscan/{sample}/{sample}.{prefix}.samtools_mpileup.call.tsv"
    threads: config["threads"]["samtools_mpileup"]["cpus"]
    run:
        shell("{samtools} mpileup -B -f {input.ref} -o {output} {input.bam} "
              "2>{log} 1>{log} ")

rule varscan_call_snp_indel:
    input:
         rules.samtools_mpileup.output
    output:
          vcf=config["dir_variants"] + "varscan/varscan_details/{sample}/{sample}.{prefix}.varscan.raw.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["varscan_call_snp_indel"]["cpus"]
    log:
       config["dir_logs"] + "varscan/{sample}/{sample}.{prefix}.call.log"
    benchmark:
             config["dir_logs"] + "varscan/{sample}/{sample}.{prefix}.call.log"

    run:
        shell("{varscan} mpileup2cns {input} {params.extra} --min-coverage 6 --min-reads2 3 "
              "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {wildcards.sample} | "
              "{bcftools} view -Oz -o {output.vcf} 1>{log} 2>{log} ")

rule varscan_combile_vcf:
    input:
         vcfs=expand(config["dir_variants"] + "varscan/varscan_details/{sample}/{sample}.{{prefix}}.varscan.raw.vcf.gz",
                     sample=config["samples"]),
         vcfs_index=expand(
             config["dir_variants"] + "varscan/varscan_details/{sample}/{sample}.{{prefix}}.varscan.raw.vcf.gz.tbi",
             sample=config["samples"])
    output:
          config["dir_variants"] + "varscan/" + config["project"] + ".{prefix}.varscan.raw.vcf.gz"
    threads: config["threads"]["varscan_combile_vcf"]["cpus"]
    log:
       config["dir_logs"] + "varscan/{prefix}.combine.log"

    benchmark:
             config["dir_logs"] + "varscan/{prefix}.combine.tsv"
    run:
        shell("{bcftools} merge -Oz --threads {threads} -o {output} {input.vcfs}")
rule varscan_marker_low_quailty:
    input:
         vcf=rules.varscan_combile_vcf.output
    output:
          vcf=config["dir_variants"] + "varscan/" + config["project"] + ".{prefix}.varscan.marked.vcf.gz"
    threads: config["threads"]["varscan_marker_low_quailty"]["cpus"]
    log:
       config["dir_logs"] + "varscan/" + config["project"] + ".{prefix}.varscan_marked_low_quality.log"
    benchmark:
             config["dir_logs"] + "varscan/" + config["project"] + ".{prefix}.varscan_marked_low_quality.tsv"
    run:
        shell("{bcftools} filter -m + -sLowQual_PJ "
              "-i ' DP<300 && DP>=10 ' "
              "-Oz -o {output.vcf} {input.vcf} 2>{log} 1>{log} ")
