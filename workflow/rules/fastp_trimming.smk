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

        rm -f {output.r1}.tmp {output.r2}.tmp {output.html}.tmp {output.json}.tmp

        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1}.tmp \
            --out2 {output.r2}.tmp \
            --html {output.html}.tmp \
            --json {output.json}.tmp \
            --thread {threads} \
            > {log} 2>&1

        mv {output.r1}.tmp {output.r1}
        mv {output.r2}.tmp {output.r2}
        mv {output.html}.tmp {output.html}
        mv {output.json}.tmp {output.json}
        """
