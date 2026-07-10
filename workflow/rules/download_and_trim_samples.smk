rule download_and_trim_fastqs:
    output:
        r1=temp(str(TRIMMED_DIR / "{sample}.trimmed.R1.fastq.gz")),
        r2=temp(str(TRIMMED_DIR / "{sample}.trimmed.R2.fastq.gz")),
        html=str(LOCAL_PATH / "logs" / "fastp" / "{sample}.html"),
        json=str(LOCAL_PATH / "logs" / "fastp" / "{sample}.json")
    params:
        r1=lambda wc: get_fastq_r1(wc.sample),
        r2=lambda wc: get_fastq_r2(wc.sample)
    threads: 4
    resources:
        fastp_limit=1,
    log:
        str(LOCAL_PATH / "logs" / "fastp" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {TRIM_DIR} {LOCAL_PATH}/logs/fastp {FASTQ_DIR}

        raw_r1="{FASTQ_DIR}/{wildcards.sample}.R1.fastq.gz"
        raw_r2="{FASTQ_DIR}/{wildcards.sample}.R2.fastq.gz"

        gcloud storage cp "{params.r1}" "$raw_r1" > {log} 2>&1
        gcloud storage cp "{params.r2}" "$raw_r2" >> {log} 2>&1

        fastp \
          --in1 "$raw_r1" \
          --in2 "$raw_r2" \
          --out1 {output.r1} \
          --out2 {output.r2} \
          --html {output.html} \
          --json {output.json} \
          --thread {threads} \
          >> {log} 2>&1

        rm -f "$raw_r1" "$raw_r2"
        """
