rule coverage_qc:
    input:
        cpg=str(
            METH_DIR / "{sample}.CpG.methylKit.gz"
        ),
        bam=str(
            BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam"
        ),
        bai=str(
            BAM_DIR / "{sample}.aligned.sorted.filt.bl.bam.bai"
        ),
        fastp_json=str(
            QC_DIR / "fastp" / "{sample}.fastp.json"
        )
    output:
        tsv=temp(
            str(COVERAGE_DIR / "{sample}.coverage_qc.tsv")
        )
    params:
        script=str(
            Path(REPO_PATH)
            / "workflow"
            / "scripts"
            / "coverage_qc.py"
        ),
        genome_size=3_100_000_000,
        min_mapq=30,
        min_baseq=0,
        excluded_contigs="lambda,pUC19"
    threads:
        4
    log:
        str(
            LOCAL_PATH
            / "logs"
            / "coverage_qc"
            / "{sample}.log"
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p {COVERAGE_DIR}
        mkdir -p "$(dirname {log})"

        python3 {params.script} \
            --sample "{wildcards.sample}" \
            --cpg "{input.cpg}" \
            --bam "{input.bam}" \
            --fastp-json "{input.fastp_json}" \
            --output "{output.tsv}" \
            --genome-size {params.genome_size} \
            --min-mapq {params.min_mapq} \
            --min-baseq {params.min_baseq} \
            --excluded-contigs "{params.excluded_contigs}" \
            > "{log}" 2>&1
        """
