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
