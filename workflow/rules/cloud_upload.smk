rule upload_results:
    input:
        raw_bam=str(BAM_DIR / "{sample}.aligned.sorted.bam"),
        raw_bai=str(BAM_DIR / "{sample}.aligned.sorted.bam.bai"),
        final_bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        final_bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        cpg=str(METH_DIR / "{sample}.CpG.methylKit.gz"),
        methyldackel_bgz=str(METH_DIR / "{sample}.methyldackel.txt.bgz"),
        methyldackel_tbi=str(METH_DIR / "{sample}.methyldackel.txt.bgz.tbi"),
        mbias_ot=str(MBIAS_DIR / "{sample}_OT.svg"),
        mbias_ob=str(MBIAS_DIR / "{sample}_OB.svg"),
        lambda_qc=str(SPIKEIN_DIR / "{sample}.lambda_qc.tsv"),
        coverage_qc=str(COVERAGE_DIR / "{sample}.coverage_qc.tsv"),
        summary=str(SUMMARY_DIR / "{sample}.qc_summary.tsv")
    output:
        done=str(UPLOAD_DIR / "{sample}.upload.done")
    params:
        bucket=config["meta"]["results_bucket"],
        prefix=config["meta"]["results_prefix"]
    log:
        str(LOCAL_PATH / "logs" / "upload_results" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log}) $(dirname {output.done})

        DEST="gs://{params.bucket}/{params.prefix}/{wildcards.sample}"

        gcloud storage cp {input.raw_bam} $DEST/ > {log} 2>&1
        gcloud storage cp {input.raw_bai} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.final_bam} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.final_bai} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.cpg} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.methyldackel_bgz} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.methyldackel_tbi} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.mbias_ot} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.mbias_ob} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.lambda_qc} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.coverage_qc} $DEST/ >> {log} 2>&1
        gcloud storage cp {input.summary} $DEST/ >> {log} 2>&1

        touch {output.done}
        """
