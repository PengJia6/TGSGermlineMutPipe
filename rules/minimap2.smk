def get_one_sample_fq(wildcards):
    return config["samples_info"][wildcards.cohort][wildcards.sample]["subsamples"][wildcards.subsample]


rule minimap2:
    input:
         fq=get_one_sample_fq,
         ref=get_ref
    output:
          config[
              "dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.minimap2.sorted.addrg.bam"
    threads:  config["threads"]["minimap2"]["cpus"]
    log: config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.minimap2.merge.logs"
    benchmark:
             config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.minimap2.merge.tsv"
    run:
        shell("{minimap2}  -a -H -k19 -O 5,56 -E 4,1 -A 2 -B 5 -z 400,50 -r 2000 -g 5000 --eqx --MD -Y "
              "-R '@RG\\tID:{wildcards.subsample}\\tSM:{wildcards.sample}' "
              "{input.ref} {input.fq} |{samtools} view -Shb -@ {threads} | "
              "{samtools} sort -@ {threads} -m 2G -T {output}_tmp -o {output} -O BAM ")
        #
        # shell("{minimap2} -x map-hifi -a  --eqx "
        #       "-R '@RG\\tID:{wildcards.subsample}\\tSM:{wildcards.sample}' "
        #       "{input.ref} {input.fq} |{samtools} view -Shb -@ {threads} | "
        #       "{samtools} sort -@ {threads} -m 2G -T {output}_tmp -o {output} -O BAM ")


def get_merge_bams(wildcards):
    subsamples = list(config["samples_info"][wildcards.cohort][wildcards.sample]["subsamples"].keys())

    return expand(
        config[
            "dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}_{subsample}.{ref_name}.{aligner}.sorted.addrg.bam",
        sample=wildcards.sample, subsample=subsamples, aligner=wildcards.aligner, cohort=wildcards.cohort,
        ref_name=wildcards.ref_name)


rule sample_bam_merge:
    input:
         get_merge_bams
    output:
          bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.bam",
    log:
       config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.merge.logs"
    benchmark:
             config["dir_logs"] + "align/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.merge.tsv"

    threads: config["threads"]["sample_bam_merge"]["cpus"]
    run:
        if len(input) < 2:
            shell("cp {input} {output.bam} ")
            shell("sleep 1")
            shell("touch -h {output.bam} 1>>{log} 2>>{log}")
            shell("echo only one bam, make soft link 1>>{log} 2>>{log}")
        else:
            shell("{samtools} merge -@ {threads} {output.bam} {input} 2>{log} 1>{log}")
