# ======================================================================================================================
# Project: Project_Human_iharbor
# Script : chineseQuartetphased.smk TODO check 
# Author : Peng Jia
# Date   : 2021.09.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================

path_vcf = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/01_variants/filtering_variants/chinese_quartet.GRCh38.CCS.pbmm2.dv.pass.vcf.gz"

# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO
rule read_phased:
    input:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.bam",
        bai=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.bam.bai",
        ref=get_ref,
        vcf=path_vcf
    output:
        bam=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.phased.bam",
        list=config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.phased.list",
    log: config["dir_aligned_reads"] + "{cohort}/{sample}/{cohort}.{sample}.{ref_name}.{tech}.{aligner}.phased.bam",
    threads: 8
    params:
        extra="",
    run:
        shell("{whatshap} haplotag -o {output.bam} -r {input.ref} --sample {wildcards.sample} "
              "--output-haplotag-list  {output.list} "
              " {input.vcf} {input.bam} 2>{log} 1>{log}")