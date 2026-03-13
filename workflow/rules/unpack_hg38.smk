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
