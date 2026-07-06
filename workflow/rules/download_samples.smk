rule download_fastqs:
    output:
        r1=temp(str(FASTQ_DIR / "{sample}.R1.fastq.gz")),
        r2=temp(str(FASTQ_DIR / "{sample}.R2.fastq.gz"))
    params:
        r1=lambda wc: get_fastq_r1(wc.sample),
        r2=lambda wc: get_fastq_r2(wc.sample)
    log:
        str(LOCAL_PATH / "logs" / "download_fastqs" / "{sample}.log")
    resources:
        upload_limit=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {FASTQ_DIR} {LOCAL_PATH}/logs/download_fastqs

        gcloud storage cp "{params.r1}" "{output.r1}" > {log} 2>&1
        gcloud storage cp "{params.r2}" "{output.r2}" >> {log} 2>&1
        """
