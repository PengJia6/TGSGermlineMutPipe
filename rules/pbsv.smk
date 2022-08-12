rule pbSV_discovery:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.bam.bai",
        ref=get_ref,
    # sindex=config["ref"]["fasta"] + ".fai"
    output:
        sig=config[
                "dir_variants"] + "{cohort}/{sample}/pbsv/contigs/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.pbsv.svsig.gz"
    params:
        extra="",
    log:
        config["dir_logs"] + "pbsv/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.pbsv.log"
    benchmark:
        config["dir_logs"] + "pbsv/{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.pbsv.tsv"
    threads: config["threads"]["pbSV_discovery"]["cpus"]
    run:
        output_pre = "/".join(str(output).split("/")[:-1])
        shell("mkdir -p {output_pre}")
        shell("{pbsv} discover --region {wildcards.contig} {input.bam} {output.sig} 2>{log} 1>{log} ")


def get_pbsv_contig_sig(wildcards):
    sigs = expand(config[
                      "dir_variants"] + "{cohort}/{sample}/pbsv/contigs/{cohort}.{sample}.{ref_name}.{suffix}.{contig}.pbsv.svsig.gz",
        cohort=wildcards.cohort,ref_name=wildcards.ref_name,suffix=wildcards.suffix,sample=wildcards.sample,
        contig=config["ref"][wildcards.ref_name]["avaliable"]
                  )
    return sigs


rule pbsv_call:
    input:
        sigs=get_pbsv_contig_sig,
        ref=get_ref
    output:
        vcf=config["dir_variants"] + "{cohort}/{sample}/pbsv/{cohort}.{sample}.{ref_name}.{suffix}.pbsv.raw.vcf.gz"
    threads: config["threads"]["pbsv_call"]["cpus"]
    run:
        vcf = str(output.vcf)[:-3]
        shell("{pbsv} call -t DEL,INS,INV,DUP,BND,CNV --ccs -j {threads} {input.ref} {input.sigs} {vcf}")
        shell("{bcftools} view -Oz -o {output.vcf} {vcf}")
