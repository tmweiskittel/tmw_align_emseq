rule sort_bam:
    input:
        bam=str(BAM_DIR / "{sample}.bam")
    output:
        bam=str(BAM_DIR / "{sample}.sorted.bam")
    threads: 8
    log:
        "logs/sort/{sample}.log"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam} \
        > {log} 2>&1
        """
rule filter_bam:
    input:
        bam=str(BAM_DIR / "{sample}.sorted.bam")
    output:
        bam=str(BAM_DIR / "{sample}.sorted.filt.bam")
    threads: 8
    log:
        "logs/filter/{sample}.filter.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        samtools view -@ {threads} -u -f 2 -q 30 -F 3840 {input.bam} \
          | samtools sort -@ {threads} -o {output.bam} - \
          > {log} 2>&1
        """
        
rule index_filtered_bam:
    input:
        bam=str(BAM_DIR / "{sample}.sorted.filt.bam")
    output:
        bai=str(BAM_DIR / "{sample}.sorted.filt.bam.bai")
    log:
        "logs/index/{sample}.sorted.filt.index.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bai}) $(dirname {log})
        samtools index {input.bam} > {log} 2>&1
        """

  rule blacklist_filter_bam:
    input:
        bam=str(BAM_DIR / "{sample}.sorted.filt.bam"),
        bai=str(BAM_DIR / "{sample}.sorted.filt.bam.bai"),
        blacklist=config["ref"]["blacklist"]
    output:
        bam=str(BAM_DIR / "{sample}.sorted.filt.bl.bam")
    threads: 8
    log:
        "logs/blacklist/{sample}.blacklist.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        bedtools intersect -abam {input.bam} -b {input.blacklist} -v \
          | samtools sort -@ {threads} -o {output.bam} - \
          > {log} 2>&1
        """

rule index_blacklist_filtered_bam:
    input:
        bam=str(BAM_DIR / "{sample}.sorted.filt.bl.bam")
    output:
        bai=str(BAM_DIR / "{sample}.sorted.filt.bl.bam.bai")
    log:
        "logs/index/{sample}.sorted.filt.bl.index.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bai}) $(dirname {log})
        samtools index {input.bam} > {log} 2>&1
        """
