rule gatk_hc_call:
    input:
         bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam",
         bai=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam.bai",
         ref=config["ref"]["fasta"],
    output:
          gvcf=config["dir_variants"] + "gatk/gatk_details/{sample}/{sample}.{prefix}.{contig}.gvcf.gz"
    log:
       config["dir_logs"] + "gatk/{sample}/{sample}.{prefix}.{contig}.call.log"
    threads: config["threads"]["gatk_hc_call"]["cpus"]
    benchmark:
             config["dir_logs"] + "gatk/{sample}/{sample}.{prefix}.{contig}.call.tsv"
    params:
          extra="",
          java_options="",
          regions="",
          dbsnp=[],
    run:
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt} HaplotypeCaller {params.extra} --minimum-mapping-quality 8  "
              " -R {input.ref} -ERC GVCF -L {wildcards.contig} -I {input.bam} -O {output.gvcf}"
              " 2>{log} 1>{log}")

rule gatk_combine_gvcf:
    input:
         ref=config["ref"]["fasta"],
         gvcfs=expand(config["dir_variants"] + "gatk/gatk_details/{sample}/{sample}.{{prefix}}.{{contig}}.gvcf.gz",
                      sample=config["samples"])
    output:
          config["dir_variants"] + "gatk/gatk_details/contigs/" + config["project"] + ".{prefix}.{contig}.gvcf.gz"
    params:
          extra="",
          java_options=""
    threads: config["threads"]["gatk_combine_gvcf"]["cpus"]
    log:
       config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.{contig}.combined.log"
    benchmark:
             config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.{contig}.combined.tsv"

    run:
        inputs = " ".join([("-V " + f) for f in input.gvcfs])
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt} CombineGVCFs {params.extra} "
              " {inputs} -R {input.ref} -O {output} 2>{log} 1>{log}")

rule gatk_genotype:
    input:
         ref=config["ref"]["fasta"],
         gvcf=rules.gatk_combine_gvcf.output
    output:
          vcf=config["dir_variants"] + "gatk/gatk_details/contigs/" + config["project"] + ".{prefix}.{contig}.vcf.gz"
    params:
          extra="",
          java_options=""
    threads: config["threads"]["gatk_genotype"]["cpus"]
    log:
       config["dir_logs"] + "gatk/" + config["project"] + ".genotype.{prefix}.{contig}.log"
    benchmark:
             config["dir_logs"] + "gatk/" + config["project"] + ".genotype.{prefix}.{contig}.tsv"
    run:
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt}  GenotypeGVCFs {params.extra} "
              " -R {input.ref} -V {input.gvcf} -O {output.vcf} 2>{log} 1>{log}")

rule gatk_merge_contig_vcf:
    input:
         expand(config["dir_variants"] + "gatk/gatk_details/contigs/" +
                config["project"] + ".{{prefix}}.{contig}.vcf.gz",
                contig=config["contigs"])
    output:
          config["dir_variants"] + "gatk/" + config["project"] + ".{prefix}.gatk.raw.vcf.gz"
    params:
          extra="",
          java_options=""
    threads: config["threads"]["gatk_merge_contig_vcf"]["cpus"]
    log:
       config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.gatk_merge_contig_vcf.log"
    benchmark:
             config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.gatk_merge_contig_vcf.tsv"
    run:
        # inputs = " ".join([("-V " + f) for f in input.gvcfs])
        inputs = " ".join(["INPUT={}".format(f) for f in input])
        shell("{picard} MergeVcfs {params.extra} "
              " {inputs} OUTPUT={output} 2>{log} 1>{log}")

rule gatk_marker_low_quailty:
    input:
         vcf=rules.gatk_merge_contig_vcf.output
    output:
          vcf=config["dir_variants"] + "gatk/" + config["project"] + ".{prefix}.gatk.marked.vcf.gz"
    threads: config["threads"]["gatk_marker_low_quailty"]["cpus"]
    log:
       config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.gatk_marked_low_quality.log"
    benchmark:
             config["dir_logs"] + "gatk/" + config["project"] + ".{prefix}.gatk_marked_low_quality.tsv"
    run:
        shell("{bcftools} filter -m + -sLowQual_PJ --threads {threads} "
              "-i 'QD>2.0 && FS<60.0 && ReadPosRankSum>-8.0 && FORMAT/DP<300 && FORMAT/DP>=10' "
              "-Oz -o {output.vcf} {input.vcf}  2>{log} 1>{log}")

# # inputs = " ".join([("-V " + f) for f in input.gvcfs])
# inputs = " ".join(["INPUT={}".format(f) for f in input])
# shell("{path_picard}picard MergeVcfs "
#       " {inputs} OUTPUT={output} 2>{log} 1>{log}")
