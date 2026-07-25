rule coverage_qc:
    input:
        cpg=str(
            METH_DIR / "{sample}.CpG.methylKit.gz"
        ),
        bam=str(
            ALIGN_DIR / "{sample}.sorted.bam"
        ),
        bai=str(
            ALIGN_DIR / "{sample}.sorted.bam.bai"
        ),
        fastp_json=str(
            FASTP_DIR / "{sample}.fastp.json"
        )
    output:
        tsv=temp(
            str(
                COVERAGE_DIR
                / "{sample}.coverage_qc.tsv"
            )
        )
    params:
        genome_size=3_100_000_000,
        min_mapq=20,
        min_baseq=0,
        excluded_contigs=[
            "lambda",
            "pUC19",
        ]
    log:
        str(
            LOCAL_PATH
            / "logs"
            / "coverage_qc"
            / "{sample}.log"
        )
    script:
        "workflow/scripts/coverage_qc.py"
