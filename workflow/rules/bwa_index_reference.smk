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
