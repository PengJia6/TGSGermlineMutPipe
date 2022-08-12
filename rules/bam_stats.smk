rule mosdepth_stats:
    input:
         bam=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam",
         bai=config["dir_aligned_reads"] + "{sample}/{sample}.{prefix}.bam.bai",

    output:
          config["dir_reports"] + "aligned_reads/{sample}/{sample}.{prefix}.mosdepth.summary.txt",

    threads: config["threads"]["mosdepth"]["cpus"]
    run:
        output_prefix = str(output).rstrip(".summary.txt")
        shell("{mosdepth} -t{threads} {output} {input.bam}")
