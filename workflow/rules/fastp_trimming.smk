rule fastp_trim:
    input:
        r1=str(FASTQ_DIR / "{sample}.R1.fastq.gz"),
        r2=str(FASTQ_DIR / "{sample}.R2.fastq.gz")
    output:
        r1=temp(str(TRIMMED_DIR / "{sample}.trimmed.R1.fastq.gz")),
        r2=temp(str(TRIMMED_DIR / "{sample}.trimmed.R2.fastq.gz")),
        html=temp(str(QC_DIR / "fastp" / "{sample}.fastp.html")),
        json=temp(str(QC_DIR / "fastp" / "{sample}.fastp.json"))
    conda:
        "../envs/fastp.yaml"
    threads: 2
    log:
        str(LOCAL_PATH / "logs" / "fastp" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {TRIMMED_DIR} {QC_DIR}/fastp {LOCAL_PATH}/logs/fastp

        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            > {log} 2>&1
        """
