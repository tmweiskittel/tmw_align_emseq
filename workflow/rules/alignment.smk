rule bwameth_align:
    input:
        r1=str(TRIMMED_DIR / "{sample}.trimmed.R1.fastq.gz"),
        r2=str(TRIMMED_DIR / "{sample}.trimmed.R2.fastq.gz"),
        ref=str(BWA_FA),
        amb=str(BWA_FA) + ".amb",
        ann=str(BWA_FA) + ".ann",
        bwt=str(BWA_FA) + ".bwt",
        pac=str(BWA_FA) + ".pac",
        sa=str(BWA_FA) + ".sa"
    output:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.bam"),
        bai=str(BAM_DIR / "{sample}.aligned.sorted.bam.bai")
    conda:
        "../envs/bwameth.yaml"    
    params:
        bwameth_threads=12,
        sort_threads=4
    threads: 16
    log:
        str(LOCAL_PATH / "logs" / "bwameth_align" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {BAM_DIR} {LOCAL_PATH}/logs/bwameth_align

        bwameth.py \
            --threads {params.bwameth_threads} \
            --reference {input.ref} \
            {input.r1} {input.r2} \
            2> {log} \
        | samtools sort \
            -@ {params.sort_threads} \
            -o {output.bam} \
            - 2>> {log}

        samtools index {output.bam} {output.bai} >> {log} 2>&1
        """
