rule vcf_normalization:
    input:
         vcf="{prefix}.vcf.gz",
         vcf_index="{prefix}.vcf.gz"
    output:
          "{prefix}.norm.vcf.gz"
    threads: config["threads"]["vcf_normalization"]["cpus"]
    run:
        shell("{bcftools} norm -m+  --threads {threads} -o {output} -Oz {input.vcf}")
