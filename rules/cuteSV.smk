rule cuteSV:
    input:
         bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
         bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
         ref=get_ref,
         # sindex=config["ref"]["fasta"] + ".fai"
    output:
          vcf=config[
                  "dir_variants"] + "{cohort}/{sample}/cuteSV/{cohort}.{sample}.{ref_name}.{suffix}.cuteSV.raw.vcf.gz"
    params:
          extra="",
    log:
       config["dir_logs"] + "cuteSV/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.cuteSV.log"
    benchmark:
             config["dir_logs"] + "cuteSV/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.cuteSV.tsv"
    threads: config["threads"]["cuteSV"]["cpus"]
    run:
        workdir = str(output.vcf).rstrip(".vcf.gz") + "_tmp"
        vcf = str(output.vcf).rstrip(".gz")
        shell("mkdir -p {workdir}")
        shell("{cuteSV} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 "
              "-s 2 --genotype "
              "--diff_ratio_merging_DEL	0.5 --max_cluster_bias_DEL	1000 "
              "--threads {threads} --sample	{wildcards.sample} --min_size 20 "
              "{input.bam} {input.ref} {vcf} {workdir} 2>{log} 1>{log} ")
        shell("{bcftools} view -Oz -o {output.vcf} {vcf}")
