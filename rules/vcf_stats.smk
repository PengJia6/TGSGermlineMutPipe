# ======================================================================================================================
# Project: Project_Human
# Script : vcf_stats.smk TODO check 
# Author : Peng Jia
# Date   : 2021.04.28
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================

# localrules: a
# ruleorder: a > b

# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO
rule vcf_qc:
    input:
         vcf="{prefix}.vcf.gz",
         vcf_index="{prefix}.vcf.gz.tbi",
         ref=get_ref
    output:
          vcf="{prefix}.VariantQC"

    params:
          extra="",
    run:
        shell("java -jar {discvrseq} -O {output.vcf} -R {input.ref} -V{input.vcf}")
        # shell("command 2")
