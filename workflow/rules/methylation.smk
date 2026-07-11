rule methyldackel_extract:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        ref=str(BWA_FA)
    output:
        cpg=temp(str(METH_DIR / "{sample}.CpG.methylKit"))
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
        mkdir -p $(dirname {output.cpg}) $(dirname {log})

        MethylDackel extract \
            -@ {threads} \
            --methylKit \
            --minDepth 5 \
            --maxVariantFrac 0.5 \
            -o {params.prefix} \
            {input.ref} \
            {input.bam} \
            > {log} 2>&1
        mv {params.prefix}_CpG.methylKit {output.cpg}
        """

rule methyldackel_mbias:
    input:
        bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        ref=str(BWA_FA)
    output:
        ot=temp(str(MBIAS_DIR / "{sample}_OT.svg")),
        ob=temp(str(MBIAS_DIR / "{sample}_OB.svg"))
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

rule make_single_methylkit_methyldackel_obj:
    input:
        cpg=str(METH_DIR / "{sample}.CpG.methylKit")
    output:
        bgz=temp(str(METH_DIR / "{sample}.methyldackel.txt.bgz")),
        tbi=temp(str(METH_DIR / "{sample}.methyldackel.txt.bgz.tbi")),
        cpg_gz=temp(str(METH_DIR / "{sample}.CpG.methylKit.gz"))
    params:
        Rscript=f"{REPO_PATH}/workflow/scripts/make_single_amp_methylkit_obj.R",
        mincov=config.get("emseq_mincov", 5),
        build=config["meta"]["ref_name"],
        treatment=1
    threads: 1
    log:
        str(LOCAL_PATH / "logs" / "make_single_methylkit_methyldackel_obj" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bgz}) $(dirname {log})

       conda run -n methylkit Rscript {params.Rscript} \
            --amp_file {input.cpg} \
            --library_id {wildcards.sample}.methyldackel \
            --mincov {params.mincov} \
            --out_dir $(dirname {output.bgz}) \
            --treatment {params.treatment} \
            --build {params.build} \
    > {log} 2>&1
        gzip -c {input.cpg} > {output.cpg_gz}
        """
