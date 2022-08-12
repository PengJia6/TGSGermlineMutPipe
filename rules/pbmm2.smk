rule pbmm2_align:
    input:
         fq=config["dir_raw_data"] + "{sample}.CCS.fofn",  #rules.common
         ref=get_ref
    output:
          # bam=dir_aligned_reads + "{sample}/{sample}.GRCh38.CCS.ngmlr.bam",
          bam=config["dir_aligned_reads"] + "{sample}/{sample}." + config["ref"]["name"] + ".CCS.pbmm2.bam",
          # bam=config["dir_aligned_reads"] + "{sample}/{sample}_{subsample}.pbmm2.sorted.addrg.bam"
    log:
       config["dir_logs"] + "align/{sample}.GRCh38.CCS.pbmm2.logs"
    benchmark:
             config["dir_logs"] + "align/{sample}.GRCh38.CCS.pbmm2.tsv"
    threads: config["threads"]["pbmm2_align"]["cpus"]
    params:
          extra="",
    run:
        shell("{pbmm2} align {input.ref} {input.fq} {output.bam} --preset CCS --unmapped "
              "--sample {wildcards.sample} "
              "--sort -m 8G --sample {wildcards.sample} -J {threads} -j {threads} 2>{log} 1>{log} ")
# shell("{pbmm2} align {input.ref} {input.fq} {output.bam} --preset CCS --unmapped "
#            "--rg '@RG\tID:{wildcards.sample}_{wildcards.subsample}\tSM:{wildcards.sample}\tLB:{wildcards.subsample}\tPL:PACBIO\tPU:unknown\tDS:CCS' "
#            "--sort -m 8G --sample {wildcards.sample} -J {threads} -j {threads} 2>{log} 1>{log} ")
