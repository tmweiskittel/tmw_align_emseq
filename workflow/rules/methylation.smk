rule methyldackel_extract:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        ref=str(BWA_FA)
    output:
        cpg=str(METH_DIR / "{sample}.aligned.sorted.filt.bl_methyldackel_CpG.methylKit")
    conda:
        "../envs/methyldackel.yaml"
    threads: 8
    log:
        str(LOCAL_PATH / "logs" / "methyldackel_extract" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {METH_DIR} $(dirname {log})

        MethylDackel extract \
            -@ {threads} \
            --methylKit \
            --mergeContext \
            --minDepth 5 \
            --maxVariantFrac 0.5 \
            -o {METH_DIR}/{wildcards.sample}.aligned.sorted.filt.bl_methyldackel \
            {input.ref} \
            {input.bam} \
            > {log} 2>&1
        """

rule methyldackel_mbias:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        ref=str(BWA_FA)
    output:
        txt=str(MBIAS_DIR / "{sample}.aligned.sorted.filt.bl_methyldackel.M-bias.txt")
    conda:
        "../envs/methyldackel.yaml"
    threads: 8
    log:
        str(LOCAL_PATH / "logs" / "methyldackel_mbias" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {MBIAS_DIR} $(dirname {log})

        MethylDackel mbias \
            -@ {threads} \
            -o {MBIAS_DIR}/{wildcards.sample}.aligned.sorted.filt.bl_methyldackel \
            {input.ref} \
            {input.bam} \
            > {log} 2>&1
        """
