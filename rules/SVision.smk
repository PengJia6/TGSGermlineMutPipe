rule SVision:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
        model="/data/home/pengjia/mysoftware/mutation/SVision/models/svision-cnn-model.ckpt.data-00000-of-00001"
    # sindex=config["ref"]["fasta"] + ".fai"
    output:
        vcf=config[
                "dir_variants"] + "{cohort}/{sample}/SVision/{cohort}.{sample}.{ref_name}.{suffix}.SVision.raw.vcf.gz"
    params:
        extra="",
    log:
        config["dir_logs"] + "SVision/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.SVision.log"
    benchmark:
        config["dir_logs"] + "SVision/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.SVision.tsv"
    threads: config["threads"]["SVision"]["cpus"]
    run:
        # workdir = str(output.vcf).rstrip(".vcf.gz") + "_tmp"
        # vcf = str(output.vcf).rstrip(".gz")
        dir = config["dir_variants"] + f"{wildcards.cohort}/{wildcards.sample}/SVision/" \
                                       f"{wildcards.cohort}.{wildcards.sample}.{wildcards.ref_name}.{wildcards.suffix}.SVision_dir"
        my_model = "/data/home/pengjia/mysoftware/mutation/SVision/models/svision-cnn-model.ckpt"
        shell("{SVision} -t {threads} -o {dir} -n {wildcards.sample} -b {input.bam} "
              "-g {input.ref} -m {my_model} --min_sv_size 30 -s10  2>{log} 1>{log}")
        shell("{bgzip} -c {dir}/{wildcards.sample}.svision.s10.vcf > {output.vcf}")
        shell("touch {output.vcf}")
#
# shell("{cuteSV} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 "
#       "-s 2 "
#       "--diff_ratio_merging_DEL	0.5 --max_cluster_bias_DEL	1000 "
#       "--threads {threads} --sample	{wildcards.sample} --min_size 20 "
#       "{input.bam} {input.ref} {vcf} {workdir} 2>{log} 1>{log} ")
# shell("{bcftools} view -Oz -o {output.vcf} {vcf}")
