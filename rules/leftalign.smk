## TODO check
rule left_align:
    input:
         bam=config["dir_aligned_reads"] + "{prefix}.bam",
         bai=config["dir_aligned_reads"] + "{prefix}.bam.bai",
         ref= config["ref"]["fasta"]
    output:
          bam=config["dir_aligned_reads"] + "{prefix}.leftAlign.bam"
    log:
       config["dir_logs"] + "leftAlign/{prefix}.logs"
    benchmark:
             config["dir_logs"] + "leftAlign/{prefix}.logs"
    threads: config["threads"]["left_align"]["cpus"]
    run:
        tmp=str(output.bam).rstrip(".bam")+"_tmp"
        shell("mkdir -p {tmp}")
        shell("{gatk} LeftAlignIndels -I {input.bam} -O {output.bam} -R {input.ref} "
              "--tmp-dir {tmp} 2>{log} 1>{log}")
