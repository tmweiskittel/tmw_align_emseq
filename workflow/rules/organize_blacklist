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
