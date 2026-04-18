rule methyldackel_extract:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        ref=str(BWA_FA)
    output:
        cpg=str(METH_DIR / "{sample}.CpG.methylKit.gz")
    params:
        prefix=str(METH_DIR / "{sample}")
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
            --minDepth 5 \
            --maxVariantFrac 0.5 \
            -o {params.prefix} \
            {input.ref} \
            {input.bam} \
            > {log} 2>&1
        gzip -f {params.prefix}_CpG.methylKit
        mv {params.prefix}_CpG.methylKit.gz {output.cpg}
        """

rule methyldackel_mbias:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        ref=str(BWA_FA)
    output:
        txt=str(MBIAS_DIR / "{sample}.M-bias.txt")
    params:
        prefix=str(MBIAS_DIR / "{sample}")
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
            {input.ref} \
            {input.bam} \
            {params.prefix} \
            > {log} 2>&1
        

        """

