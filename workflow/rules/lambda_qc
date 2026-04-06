rule lambda_spikein_qc:
    input:
        cpg=str(METH_DIR / "{sample}.CpG.methylKit.gz")
    output:
        tsv=str(SPIKEIN_DIR / "{sample}.lambda_qc.tsv")
    log:
        str(LOCAL_PATH / "logs" / "lambda_spikein_qc" / "{sample}.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SPIKEIN_DIR} $(dirname {log})

        python - <<'PY' > {log} 2>&1
import gzip
import csv

cpg_file = "{input.cpg}"
out_file = "{output.tsv}"

sites = 0
meth = 0
unmeth = 0

with gzip.open(cpg_file, "rt") as fh:
    for line in fh:
        if line.startswith("chr"):
            continue

        fields = line.strip().split("\t")
        chrom = fields[0]

        if chrom != "lambda":
            continue

        m = int(fields[5])
        u = int(fields[6])

        sites += 1
        meth += m
        unmeth += u

mean_methyl = "NA"
if meth + unmeth > 0:
    mean_methyl = meth / (meth + unmeth)

with open(out_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "lambda_cpg_sites",
        "lambda_methylated_counts",
        "lambda_unmethylated_counts",
        "lambda_mean_methylation_fraction"
    ])
    writer.writerow([sites, meth, unmeth, mean_methyl])
PY
        """
