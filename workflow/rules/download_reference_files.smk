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
