from pathlib import Path

SAMPLESHEET = Path(DATA_PATH) / SAMPLE_SHEET
MD5_FILE = Path(DATA_PATH) / MD5_SHEET

rule unpack_hg38:
    input:
        str(HG38_GZ)
    output:
        str(HG38_FA)
    log:
        "logs/unpack_hg38.log"
    shell:
        r"""
        set -euo pipefail
        gunzip -c {input} > {output} 2> {log}
        """

rule organize_blacklist:
    input:
        str(BLACKLIST_GZ)
    output:
        str(BLACKLIST_BED_GZ)
    log:
        "logs/organize_blacklist.log"
    shell:
        r"""
        set -euo pipefail
        cp {input} {output} 2> {log}
        """

rule normalize_spikein_headers:
    input:
        lambda_in=str(LAMBDA_RAW_FA),
        puc19_in=str(PUC19_RAW_FA)
    output:
        lambda_out=str(LAMBDA_RENAMED_FA),
        puc19_out=str(PUC19_RENAMED_FA)
    log:
        "logs/normalize_spikein_headers.log"
    shell:
        r"""
        set -euo pipefail

        awk 'BEGIN{{printed=0}} /^>/ {{if(!printed){{print ">lambda"; printed=1}}; next}} {{print}}' \
            {input.lambda_in} > {output.lambda_out} 2> {log}

        awk 'BEGIN{{printed=0}} /^>/ {{if(!printed){{print ">pUC19"; printed=1}}; next}} {{print}}' \
            {input.puc19_in} > {output.puc19_out} 2>> {log}
        """

rule download_reference_files:
    output:
        hg38_gz=str(HG38_GZ),
        blacklist_gz=str(BLACKLIST_GZ),
        lambda_fa=str(LAMBDA_RAW_FA),
        puc19_fa=str(PUC19_RAW_FA)
    params:
        hg38_url=config["reference"]["hg38_fasta_gz"],
        blacklist_url=config["reference"]["hg38_blacklist_bed_gz"],
        lambda_url=config["reference"]["lambda_fasta"],
        puc19_url=config["reference"]["puc19_fasta"]
    log:
        "logs/download_reference_files.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {REF_SOURCE} {REF_FASTA} {REF_BED} {REF_BWA} logs

        curl -L "{params.hg38_url}" -o {output.hg38_gz} > {log} 2>&1
        curl -L "{params.blacklist_url}" -o {output.blacklist_gz} >> {log} 2>&1
        curl -L "{params.lambda_url}" -o {output.lambda_fa} >> {log} 2>&1
        curl -L "{params.puc19_url}" -o {output.puc19_fa} >> {log} 2>&1
        """

checkpoint download_metadata:
    output:
        samples=str(LOCAL_SAMPLES),
        md5=str(LOCAL_MD5)
    params:
        samples=str(SAMPLESHEET),
        md5=str(MD5_FILE)
    log:
        str(LOCAL_PATH / "logs" / "download_metadata.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {LOCAL_PATH} {LOCAL_PATH}/logs
        gcloud storage cp gs://{params.samples} {output.samples} > {log} 2>&1
        gcloud storage cp gs://{params.md5} {output.md5} >> {log} 2>&1
        """

rule bwa_index_reference:
    input:
        str(COMBINED_FA)
    output:
        fa=str(BWA_FA),
        amb=str(BWA_FA) + ".amb",
        ann=str(BWA_FA) + ".ann",
        bwt=str(BWA_FA) + ".bwt",
        pac=str(BWA_FA) + ".pac",
        sa=str(BWA_FA) + ".sa"
    conda:
        "../envs/bwa.yaml"
    threads: 4
    log:
        "logs/bwa_index_reference.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {REF_BWA}
        cp {input} {output.fa}
        bwa index {output.fa} > {log} 2>&1
        """

rule build_combined_reference:
    input:
        hg38=str(HG38_FA),
        lambda_fa=str(LAMBDA_RENAMED_FA),
        puc19_fa=str(PUC19_RENAMED_FA)
    output:
        str(COMBINED_FA)
    log:
        "logs/build_combined_reference.log"
    shell:
        r"""
        set -euo pipefail
        cat {input.hg38} {input.lambda_fa} {input.puc19_fa} > {output} 2> {log}
        """

