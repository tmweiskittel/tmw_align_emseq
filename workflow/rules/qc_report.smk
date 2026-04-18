rule sample_qc_summary:
    input:
        fastp_json=str(QC_DIR / "fastp" / "{sample}.fastp.json"),
        raw_bam=str(BAM_DIR / "{sample}.aligned.sorted.bam"),
        raw_bai=str(BAM_DIR / "{sample}.aligned.sorted.bam.bai"),
        final_bam=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"),
        final_bai=str(BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"),
        cpg=str(METH_DIR / "{sample}.CpG.methylKit.gz"),
        lambda_qc=str(SPIKEIN_DIR / "{sample}.lambda_qc.tsv"),
        coverage_qc=str(COVERAGE_DIR / "{sample}.coverage_qc.tsv")
    output:
        tsv=temp(str(SUMMARY_DIR / "{sample}.qc_summary.tsv"))
    log:
        str(LOCAL_PATH / "logs" / "sample_qc_summary" / "{sample}.log")
    conda:
        "../envs/qc_py.yaml"
    script:
        "../scripts/sample_qc_summary.py"
