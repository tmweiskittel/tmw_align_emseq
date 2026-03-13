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
