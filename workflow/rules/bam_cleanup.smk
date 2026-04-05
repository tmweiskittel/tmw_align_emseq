rule filter_bam:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.bam")
    output:
        bam=temp(str(BAM_DIR / "{sample}.aligned.sorted.filt.bam"))
    threads: 8
    log:
        str(LOCAL_PATH / "logs" / "filter_bam" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {BAM_DIR} {LOCAL_PATH}/logs/filter_bam

        samtools view -@ {threads} -u -f 2 -q 30 -F 3840 {input.bam} \
          | samtools sort -@ {threads} -o {output.bam} - \
          > {log} 2>&1
        """

rule blacklist_filter_bam:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bam"),
        blacklist=str(BLACKLIST_BED)
    output:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam")
    threads: 8
    log:
        str(LOCAL_PATH / "logs" / "blacklist_filter_bam" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {BAM_DIR} $(dirname {log})

        bedtools intersect -abam {input.bam} -b {input.blacklist} -v \
          | samtools sort -@ {threads} -o {output.bam} - \
          > {log} 2>&1
        """

rule index_blacklist_filtered_bam:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam")
    output:
        bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai")
    log:
        str(LOCAL_PATH / "logs" / "index_blacklist_filtered_bam" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {LOCAL_PATH}/logs/index_blacklist_filtered_bam

        samtools index {input.bam} {output.bai} > {log} 2>&1
        """
